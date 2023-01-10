% Author: Lukas Ma
% 
% Date: 2023-01-10
%
% Plot figures showing the non-linear behaviour of unemployment, job
% finding probability, and uncertainty across different economic regimes

%% INITIALISING
clear all; 
clear session; 
close all; clc
warning off all

% Define working directory
pwd
currentFolder = pwd;

folder = './Results'; %Specify output folder
%%
data_lib = readtable('./Data/data_quarterly_merged.csv');

%%
% Specify starting period
% Full sample 

StartDate = datetime(1976,01,01); %Start date: Earliest record of job finding probability (1976Q1)

EndDate = datetime(2022,04,01); %End date: Latest record of GDP (2022Q2)

% Trim the data
data_lib(data_lib.Date<StartDate,:)=[]; %Remove observations before the start date
data_lib(data_lib.Date>EndDate,:)=[]; %Remove observations after the end date

%%

% Define an array of recession (monthly)
CAN_peak_m = ['1981-06-01';'1990-03-01';'2008-10-01';'2020-02-01'];

CAN_trough_m = ['1982-10-01';'1992-05-01';'2009-05-01';'2020-04-01'];

% Percentage of time the economy is in recession according to C.D. Howe
% institute Data_BCC_Aug2021.xlsx

% Define an array of recession (quarter)

CAN_peak_q = ['1981-07-01';'1990-04-01';'2008-10-01';'2020-01-01']; 

CAN_trough_q = ['1983-01-01';'1992-07-01';'2009-07-01';'2020-04-01']; %Adjusted to include the quarter in which trough occured. E.g. through occured in 1982Q4 -> Adjust by one quarter to 1983Q1


CAN_recession_q = [datetime(CAN_peak_q), datetime(CAN_trough_q)];

CAN_recession_m = [datetime(CAN_peak_m), datetime(CAN_trough_m)];

%% Insepction plot

figure(1);

%Unemployment 
subplot(1,3,1)
plot(data_lib.Date,data_lib.UNEMP_CAN, '-',LineWidth=1.5) %using Canadian business cycle dating
recessionplot('recessions',CAN_recession_q)
title('Unemployment (in percentage)')

%Job finding probability
subplot(1,3,2)
plot(data_lib.Date,data_lib.mean_PROB_JOB_FINDING,'-',LineWidth=1.5)
recessionplot('recessions',CAN_recession_q)
title('Job finding probability')

%Uncertainty 
subplot(1,3,3)
plot(data_lib.Date, data_lib.CAN_unc_h1,'-',LineWidth=1.5)
recessionplot('recessions',CAN_recession_q)
title('Uncertainty (one quarter ahead)')

%%
% Specify starting period
% Truncated sample

StartDate = datetime(1982,01,01); %Start date: Earliest record of the uncertainty index (1982Q1)

EndDate = datetime(2020,01,01); %End date: Pre-Covid (2020Q1)

% Trim the data
data_lib(data_lib.Date<StartDate,:)=[]; %Remove observations before the start date
data_lib(data_lib.Date>EndDate,:)=[]; %Remove observations after the end date

%% Insepction plot

figure(2);

%Unemployment 
subplot(1,3,1)
plot(data_lib.Date,data_lib.UNEMP_CAN, '-',LineWidth=1.5) %using Canadian business cycle dating
recessionplot('recessions',CAN_recession_q)
title('Unemployment (in percentage)')

%Job finding probability
subplot(1,3,2)
plot(data_lib.Date,data_lib.mean_PROB_JOB_FINDING,'-',LineWidth=1.5)
recessionplot('recessions',CAN_recession_q)
title('Job finding probability')

%Uncertainty 
subplot(1,3,3)
plot(data_lib.Date, data_lib.CAN_unc_h1,'-',LineWidth=1.5)
recessionplot('recessions',CAN_recession_q)
title('Uncertainty (one quarter ahead)')

