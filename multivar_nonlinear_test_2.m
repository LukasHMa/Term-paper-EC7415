% Author: Lukas Ma
% 
% Date: 2023-01-09
%
% Script performing the multivariate non-linearity test proposed by Teräsvirta & Yang (2014) for baseline
% specification in Eq.(6)

%% PRELIMINARIES
% =======================================================================
clear all; 
clear session; 
close all; clc
warning off all

% Import relevant toolboxes
addpath('./mm_2_2/_tbx/var_tbx') 

%% Import data
%data_lib = readtable("./Data/data_quarterly_clean_full.csv"); % load full sample
data_lib = readtable("./Data/data_quarterly_clean_precovid.csv"); %load precovid data


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

% Compound growth
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

%%
%Specification 2 (replace unemployment with job finding probability)
vars = {'Unc_US','y','pr','CPI','Unc_CAN'}; % only for labelling, data is loaded directly
printvars = {'US Uncertainty','Real GDP growth','Job finding probability','Inflation (CPI)','CAN Uncertainty'}; 

data = table2array(data_lib(:,["mean_US_unc_h1","Real_GDP_CAN","mean_PROB_JOB_FINDING","CPI_CAN","CAN_unc_h1"])); 

%% Parameterisation

var_par.c_case = 1; %Estimate with constant (=1), Estimate with a trend (=2)
var_par.p = 1; %one lag
n = length(vars);
ident = 'chol';
shockpos = 1; %US uncertainty shock, CAN uncertainty shock = 5
shocksize = 0; %Size of shock: 1 StD
state.nonlinear = 'yes';
state.logistic = 'yes';
state.interacted = 'no';
state.statevar = 'MA7';
state.cq = 50; 
nboot = 10000; %The first 20% will be discarded
alpha = 90;

exdata = []; 
state.s = data_lib.MA_GDP_7_NORM_Growth; %transition variable

%% Multivariate Lagrange Multiplier test for non-linearity

[nonlinLMresults] = nonlin_testSTVAR(data, var_par, exdata, state, false, alpha); 

disp(nonlinLMresults)