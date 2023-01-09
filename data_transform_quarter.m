% Author: Lukas Ma
% 
% Version date: 2023-01-09
%
% Performs transformation and cleaning of the raw data

%% INITIALISING
clear all; 
clear session; 
close all; clc
warning off all

% Define working directory
pwd
currentFolder = pwd;

%% Data treatments (logarithm, FD etc.)

data_raw = readtable('./Data/data_quarterly_merged.csv');

% Subset data
data_rate = data_raw(:,["Date","UNEMP_CAN","mean_PROB_JOB_FINDING","mean_OUTFLOW_RATE","CAN_unc_h1","mean_US_unc_h1","CAN_unc_h4","mean_US_unc_h12"]); %Rate data (compute quarter-to-quarter growth rate, Log-diff for unemployment and job-finding probability)

data_level = data_raw(:,["Date","Nom_GDP_CAN","Real_GDP_CAN","CPI_CAN"]); %Level data (Log-FD transformation and quarter-on-quarter growth rate)

data_no_tr = data_raw(:,["Date","mean_EPU_CAN","CAN_unc_h2","WUICAN","mean_US_unc_h3"]); %No transformation (Some uncertainty indices, save for robustness check eventually)

%%
D1 = LagOp({1,-1},'Lags',[0,1]); %Define the lag operator A(L) = 1-1*L

% Compute compound growth (unemployment and job finding probability)
temp_data = 100*log(table2array(data_rate(:,2:3))); %Make a copy and take the log*100

temp_fd = [];
for j = 1:size(temp_data ,2)
    temp_fd(:,j) = filter(D1,temp_data(:,j)); %Generate FD transformed data
end 

% Compute discrete growth (all else)

temp_data2 = table2array(data_rate(:,4:end));

temp_fd2 = [];
for k = 1:size(temp_data2,2)
    temp_fd2 = [temp_fd2, 100*(temp_data2(2:end,k)-temp_data2(1:end-1,k))./temp_data2(1:end-1,k)];
end

temp_time = data_rate(2:end,1); %Extract time values, drop the first row

data_rate_tr = [temp_time array2table(temp_fd) array2table(temp_fd2)]; %Concatenate 

data_rate_tr.Properties.VariableNames = data_rate.Properties.VariableNames; %Use the same variable name as the original data

%%

data_level.Nom_GDP_CAN_Growth = data_level.Nom_GDP_CAN; 
data_level.Real_GDP_CAN_Growth = data_level.Real_GDP_CAN;

% Compute compound growth
temp_data = 100*log(table2array(data_level(:,2:4))); %Make a copy and take the log*100

temp_fd = [];
for j = 1:size(temp_data ,2)
    temp_fd(:,j) = filter(D1,temp_data(:,j)); %Generate FD transformed data
end 

temp_time = data_level(2:end,1); %Extract time values, drop the first row 

% Compute discrete growth
temp_data2 = table2array(data_level(:,5:end));

temp_fd2 = [];
for k = 1:size(temp_data2,2)
    temp_fd2 = [temp_fd2, 100*(temp_data2(2:end,k)-temp_data2(1:end-1,k))./temp_data2(1:end-1,k)];
end

data_level_tr = [temp_time array2table(temp_fd) array2table(temp_fd2)]; %Concatenate 

data_level_tr.Properties.VariableNames = data_level.Properties.VariableNames; %Use the same variable name as the original data

data_level_tr(any(isnan(data_level_tr.Nom_GDP_CAN),2),:) = []; %Remove NaNs


%% Compute transition variables

% Compound growth

% Calculate the seven quarter (backward) moving averages, i.e. MA(7)
data_level_tr.MA_GDP_7 = movmean(data_level_tr.Real_GDP_CAN,[7 0]);
data_level_tr.MA_GDP_7(1:7)= NaN;% Set the initial 7 quarters to NaN

% Calculate the 12 month (4 quaters) moving average
data_level_tr.MA_GDP_4 = movmean(data_level_tr.Real_GDP_CAN,[4 0]);
data_level_tr.MA_GDP_4(1:4)= NaN; % Set the initial 4 quarters to NaN

% Calculate the seven quarter (Centered) moving averages
data_level_tr.MA_cen_GDP_7 = movmean(data_level_tr.Real_GDP_CAN,7);
data_level_tr.MA_cen_GDP_7(1:4)= NaN; % Set the initial 4 quarters to NaN

% Calculate the seven quarter (backward) moving averages of DETREND GDP
% growth
data_level_tr.Real_GDP_CAN_detrend = detrend(data_level_tr.Real_GDP_CAN);
data_level_tr.MA_GDP_7_detrend = movmean(data_level_tr.Real_GDP_CAN_detrend,[7 0]);
data_level_tr.MA_GDP_7_detrend(1:7)= NaN;% Set the initial 7 quarters to NaN

%-------------------------------------------------------------------------------
% Quarter-on-Quarter growth

% Calculate the seven quarter (backward) moving averages, i.e. MA(7)
data_level_tr.MA_GDP_7_Growth = movmean(data_level_tr.Real_GDP_CAN_Growth,[7 0]);
data_level_tr.MA_GDP_7_Growth(1:7)= NaN;% Set the initial 7 quarters to NaN

% Calculate the 12 month (4 quaters) moving average
data_level_tr.MA_GDP_4_Growth = movmean(data_level_tr.Real_GDP_CAN_Growth,[4 0]);
data_level_tr.MA_GDP_4_Growth(1:4)= NaN; % Set the initial 4 quarters to NaN

% Calculate the seven quarter (Centered) moving averages
data_level_tr.MA_cen_GDP_7_Growth = movmean(data_level_tr.Real_GDP_CAN_Growth,7);
data_level_tr.MA_cen_GDP_7_Growth(1:4)= NaN; % Set the initial 4 quarters to NaN

% Calculate the seven quarter (backward) moving averages of DETREND GDP
% growth
data_level_tr.Real_GDP_CAN_detrend_Growth = detrend(data_level_tr.Real_GDP_CAN_Growth);
data_level_tr.MA_GDP_7_detrend_Growth = movmean(data_level_tr.Real_GDP_CAN_detrend_Growth,[7 0]);
data_level_tr.MA_GDP_7_detrend_Growth(1:7)= NaN;% Set the initial 7 quarters to NaN

%% Merge the transformed data

dfs = {data_rate_tr, data_level_tr data_no_tr};
df_merge_tr = dfs{1};

for i = 2:numel(dfs)
    df_merge_tr = outerjoin(df_merge_tr, dfs{i}, 'MergeKeys',true); 
end

%% Some housekeeping & export

StartDate = datetime(1982,04,01); %Start date: Earliest record of the uncertainty index (1982Q1)

EndDate = datetime(2022,04,01); %End date: Latest record of GDP (2022Q2)

% Trim the data
df_merge_tr(df_merge_tr.Date<StartDate,:)=[]; %Remove observations before the start date
df_merge_tr(df_merge_tr.Date>EndDate,:)=[]; %Remove observations after the end date

df_merge_tr(any(isnan(df_merge_tr.CPI_CAN),2),:) = []; %Remove NaNs

writetable(df_merge_tr,'./Data/data_quarterly_clean_full.csv') %Full sample

%% Pre-Covid data

EndDate_pre = datetime(2020,01,01); %End date: Pre-Covid (2020Q1)

% Trim the data
df_merge_tr(df_merge_tr.Date>EndDate_pre,:)=[]; %Remove observations after the end date

df_merge_tr(any(isnan(df_merge_tr.CPI_CAN),2),:) = []; %Remove NaNs

writetable(df_merge_tr,'./Data/data_quarterly_clean_precovid.csv') %Full sample
