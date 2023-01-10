% Author: Lukas Ma
% 
% Date: 2023-01-10
%
% This script performs the following: 
% (1) Calibration of smoothing parameter gamma 
% (2) Selection of lags
% (3) Estimation of linear VAR
% (4) Estimation of Smooth-transitioning VAR, using an alternative ordering that corresponds to Eq.(5)
%% PRELIMINARIES
% =======================================================================
clear all; 
clear session; 
close all; clc
warning off all

% Import relevant toolboxes
addpath('./mm_2_2/_tbx/var_tbx') 
addpath('./mm_2_2/_tbx/stvar_tbx') 
addpath('./mm_2_2/_tbx/supportfct') 
addpath('./mm_2_2/_tbx/xtratools')

folder = './Results'; %Specify output folder
%% Import data
%data_lib = readtable("./Data/data_quarterly_clean_full.csv"); % load full sample
data_lib = readtable("./Data/data_quarterly_clean_precovid.csv"); %load precovid data

%% Preparation

% Percentage of time the economy is in recession according to C.D. Howe
% institute Data_BCC_Aug2021.xlsx

% Define an array of recession (quarter)

CAN_peak_q = ['1982-04-01';'1990-04-01';'2008-10-01';'2020-01-01']; %Use the 1982Q2 as starting period

CAN_trough_q = ['1983-01-01';'1992-07-01';'2009-07-01';'2020-04-01']; %Adjusted to include the quarter in which trough occured. E.g. through occured in 1982Q4 -> Adjust by one quarter to 1983Q1

CAN_recession_q = [datetime(CAN_peak_q), datetime(CAN_trough_q)];

quarters = round(days(CAN_recession_q(1:4,2) - CAN_recession_q(1:4,1)) ./ 91); % Calculate the number of quarters in recession, plus 3 to include the quarter in which through occured. 

share_recession = (sum(quarters)-1)/size(data_lib,1); % Calculate the share of the time the economy is in recession (Pre-covid)

%share_recession = sum(quarters)/size(data_lib,1); % Calculate the share of the time the economy is in recession (full-sample)



%% Calibration of the transition variable

% Normalisation (according to Züllig)
z_var = [data_lib.MA_GDP_4,data_lib.MA_GDP_7,data_lib.MA_cen_GDP_7,data_lib.MA_GDP_7_detrend,...
    data_lib.MA_GDP_4_Growth,data_lib.MA_GDP_7_Growth,data_lib.MA_cen_GDP_7_Growth,data_lib.MA_GDP_7_detrend_Growth];

z_var_NORM = [];

for k=1:size(z_var,2)
    z_orig = z_var(:,k);
    z_orig_cons = prctile(z_orig, 50); %demean using median
    z_orig_std = std(z_orig(~isnan(z_orig)));
    z_var_NORM(:,k) = (z_orig-z_orig_cons) / z_orig_std ; 
end

% Compounded growth
data_lib.MA_GDP_4_NORM = z_var_NORM(:,1);
data_lib.MA_GDP_7_NORM = z_var_NORM(:,2);
data_lib.MA_GDP_7_cen_NORM = z_var_NORM(:,3);
data_lib.MA_GDP_7_detrend_NORM = z_var_NORM(:,4);

% Quarter-on-Quarter growth
data_lib.MA_GDP_4_NORM_Growth = z_var_NORM(:,5);
data_lib.MA_GDP_7_NORM_Growth = z_var_NORM(:,6);
data_lib.MA_GDP_7_cen_NORM_Growth = z_var_NORM(:,7);
data_lib.MA_GDP_7_detrend_NORM_Growth = z_var_NORM(:,8);

state.s_multi = z_var_NORM; 

clear z_var_NORM z_var


% Plot the normalised transition variable
figure(1)
plot(data_lib.Date, data_lib.MA_GDP_4_NORM_Growth,'-',LineWidth=2)
hold on
plot(data_lib.Date, data_lib.MA_GDP_7_NORM_Growth,'--',LineWidth=2)
hold on
plot(data_lib.Date, data_lib.MA_GDP_7_cen_NORM_Growth,':',LineWidth=2)
hold on
plot(data_lib.Date, data_lib.MA_GDP_7_detrend_NORM_Growth,'-.',LineWidth=2)
hold off

%% Use a while-loop to calibrate the smoothness parameter gamma

maxiter = 100000;
tol = 1e-5;

count = 1; %Counter
gamma_c = 0; %Initial value of gamma is set to zero
while (count < maxiter)

g0 = gamma_c; %Smoothness parameter
z = data_lib.MA_GDP_7_NORM_Growth; %Transition variable
%z = data_lib.MA_GDP_4_NORM_Growth; % Alternative transition variable
%z = data_lib.MA_GDP_7_NORM; % Alternative transition variable

F_z = exp(-g0*z)./(1+exp(-g0*z)); % Prob(Recession = 1) is modelled using logistic regression

indices = find(F_z >= 1-share_recession);

share_recession_est = size(indices,1)/size(F_z,1);

diff = abs(share_recession_est-share_recession); 
if diff < tol %termination condition
    break
else
    gamma_c = g0 + 0.0001;
end

count = count + 1; %counter
end

% Plot the transition function using the calibrated gamma
figure(4);
subplot(1,2,1)
plot(data_lib.Date, data_lib.MA_GDP_7, '-k', LineWidth=2)
recessionplot('recessions',CAN_recession_q)
title('Transition variable: MA(7) of real GDP growth rate')

subplot(1,2,2)
plot(data_lib.Date, F_z, LineWidth=2)
recessionplot('recessions',CAN_recession_q)
title('Transition function')

%% Robustness checks 
%Specification 3a (Canadian uncertainty is ordered second)
vars = {'Unc_US','Unc_CAN','y','u','CPI'}; % only for labelling, data is loaded directly
printvars = {'US Uncertainty','CAN Uncertainty','Real GDP growth','Unemployment rate','Inflation (CPI)'}; 

data = table2array(data_lib(:,["mean_US_unc_h1","CAN_unc_h1","Real_GDP_CAN","UNEMP_CAN","CPI_CAN"])); 

%% Parameterization 

h = 20; % Number of horizons (to generate IRFs) Set to 5 years
c_case = 1; %Estimate with constant (=1), Estimate with a trend (=2)

% VAR parameters
n = length(vars);
ident = 'chol';
%shockpos = 1; %US uncertainty shock = 1 (Baseline ordering)
shockpos = 2;  %CAN uncertainty shock = 2 (Alternative ordering 3ab)

shocksize = 0; %Size of shock: 1 StD

% Transition variable settings
state.nonlinear = 'yes';
state.logistic = 'yes';
state.interacted = 'no';
state.statevar = 'MA7';
%state.statevar = 'MA4';

state.gamma = gamma_c;
state.cq = 50; 

state.s = z; %transition variable
state.Fs = F_z; %transition probability

%data(any(isnan(data),2),:) = []; %Remove NaNs, if EPU is used
%state.s = z(16:end); %If EPU is used
%state.Fs = F_z(16:end); %If EPU is used  

% MCMC settings
nboot = 30000; %The first 20% will be discarded
alpha = 90;
%alpha = 70; %alternative  

exdata = []; 

%% Lag selection criteria (Based on Züllig's code)
t = size(data, 1);

%pmax = round(12*((t/100)^(1/4))); %Schwart's rule of thumb 
pmax = 8; %Set pmax equal to 8 (Lags for two years)


for p=1:pmax  
  [AIC,BIC,HQC]=InfCriteria(data, p, c_case, exdata);
   AICs(p) = AIC;
   BICs(p) = BIC;
   HQCs(p) = HQC;
end

tab = array2table([(1:pmax)',AICs',BICs',HQCs']); %Display information criteria
tab.Properties.VariableNames = {'Lags','AIC','BIC','HQC'};
disp(tab)

writetable(tab,'./Results/InfC_robust_1.xlsx') 

% Display minimum
[~,pref_AIC] = min(AICs);
[~, pref_BIC] = min(BICs);
[~, pref_HQC] = min(HQCs);
tab2 = array2table([pref_AIC, pref_BIC, pref_HQC]);
tab2.Properties.VariableNames = {'AIC', 'BIC', 'HQC'};
disp('Prefered lag order:')
disp(tab2)

% Choose the prefered lag order = lowest IC score
min_IC = min([min(AICs),min(BICs),min(HQC)]);
disp(min_IC);
%%
p = pref_AIC; %Choose lag suggested by AIC, since it has lowest IC score
clear pref_AIC pref_BIC pref_HQC AIC BIC HQC AICs BICs HQCs

%% Estimate Linear VAR
VAR = estimateVAR(data, p, c_case, exdata);
VAR.C = dyn_multipliers(VAR, h); % not identified

% Bootstrap
% organization
VAR.ident = ident;
VAR.shock = zeros(n,1);
VAR.shock(shockpos) = 1; %US uncertainty shock, CAN uncertainty shock = 5

% Cholesky decomposition (already done in estimation function)
disp(VAR.S)
eps = VAR.u*inv(VAR.S);
VAR.eps = eps(:, shockpos);
clear eps
% impulse responses
VAR.IRF = zeros(h,n);
for hh=1:h
    VAR.IRF(hh,:) = (VAR.C(:,:,hh)*VAR.S*VAR.shock)';
end

[VAR.IRFbands] = bootstrapVAR(VAR, nboot, alpha, 'residual');

%Plot impulse responses and save VAR structure
figure(5)
plotirf1(VAR.IRF, VAR.IRFbands, printvars, strcat(folder,'/lin'))

%% Estimate STVAR
STVAR = estimateSTVAR(data, state, shockpos, shocksize, p, c_case, exdata, h, nboot, alpha);

% Plot impulse responses and save VAR structure
figure(6)
plotirf2(STVAR.IRF1, STVAR.IRF1bands, STVAR.IRF2, STVAR.IRF2bands, printvars, strcat(folder,'/nonlin'), {'Expansion','Recession'}, 'Northeast')

%% Plot impulse responses all together
figure(7)
plotirf3_mod(STVAR.IRF1, STVAR.IRF1bands, STVAR.IRF2, STVAR.IRF2bands, VAR.IRF, VAR.IRFbands, printvars, strcat(folder,'/lin_nonlin'),{'Expansion','Recession','Linear'}, 6)

%% Save output

%Specification 3a

%save(strcat(folder, '/out_3a_us.mat'), 'STVAR', 'VAR','printvars','folder') %US Uncertainty shock
save(strcat(folder, '/out_3a_can.mat'), 'STVAR', 'VAR','printvars','folder') %CAN Uncertainty shock