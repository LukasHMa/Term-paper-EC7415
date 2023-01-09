% Author: Lukas Ma
% 
% Date: 2023-01-09
%
% Data preparation using the raw data.

%% INITIALISING
clear all; 
clear session; 
close all; clc
warning off all

% Define working directory
pwd
currentFolder = pwd;

%% Import data: Quarterly Nominal GDP data for Canada (FRED)
df_nomGDP_CAN = readtable('./Data/NGDPSAXDCCAQ.csv', 'ReadRowNames', 0);

df_nomGDP_CAN = renamevars(df_nomGDP_CAN,["DATE", "NGDPSAXDCCAQ"],["Date","Nom_GDP_CAN"]);

df_nomGDP_CAN.Date_str = datestr(df_nomGDP_CAN.Date,'yyyy-mm-dd'); %Create a string for merge

%% Import data: Quarterly Real GDP data for Canada (Transition variable)

df_realGDP_CAN = readtable('./Data/NGDPRSAXDCCAQ.csv', 'ReadRowNames', 0);

df_realGDP_CAN = renamevars(df_realGDP_CAN,["DATE", "NGDPRSAXDCCAQ"],["Date","Real_GDP_CAN"]);

df_realGDP_CAN.Date_str = datestr(df_realGDP_CAN.Date,'yyyy-mm-dd'); %Create a string for merge

%% Import data: Quarterly CPI Index for Canada 2015=100 (FRED)

df_CPI_CAN = readtable('./Data/CANCPIALLQINMEI.csv', 'ReadRowNames', 0);

df_CPI_CAN = renamevars(df_CPI_CAN,["DATE", "CANCPIALLQINMEI"],["Date","CPI_CAN"]);

df_CPI_CAN.CPI_CAN = df_CPI_CAN.CPI_CAN(1:end);

df_CPI_CAN.Date_str = datestr(df_CPI_CAN.Date, 'yyyy-mm-dd'); 

%% Import data: World Uncertainty Index Canada
df_WUI_CAN = readtable('./Data/WUICAN.csv', 'ReadRowNames',0);

df_WUI_CAN = renamevars(df_WUI_CAN, "DATE", "Date");

df_WUI_CAN.Date_str = datestr(df_WUI_CAN.Date, 'yyyy-mm-dd');

%% Import data: Unemployment rate Canada (FRED) 
df_UNEMP_CAN = readtable('./Data/LRUNTTTTCAQ156S.csv','ReadRowNames',0);

df_UNEMP_CAN = renamevars(df_UNEMP_CAN, ["DATE","LRUNTTTTCAQ156S"], ["Date","UNEMP_CAN"]);

df_UNEMP_CAN.UNEMP_CAN = df_UNEMP_CAN.UNEMP_CAN(1:end);

df_UNEMP_CAN.Date_str = datestr(df_UNEMP_CAN.Date,'yyyy-mm-dd');


%% Import Uncertainty data (Quarterly): 
df_MacroUN_CAN = readtable('./Data/uncertainty_series_quarterly-1.xlsx', 'ReadRowNames',0); 

df_MacroUN_CAN = df_MacroUN_CAN(:,["Time","CAN_unc_h1","CAN_unc_h2","CAN_unc_h4"]);

df_MacroUN_CAN = renamevars(df_MacroUN_CAN,"Time","Date");

df_MacroUN_CAN.Date_str = datestr(df_MacroUN_CAN.Date,'yyyy-mm-dd');

%% Import data: Unemployment duration (StatCan Table: 14-10-0342-01)
df_UNEMPDUR_CAN = readtable('./Data/1410034201_databaseLoadingData.csv', 'ReadRowNames',0);

% Manipulation
df_UNEMPDUR_CAN = df_UNEMPDUR_CAN(:,["REF_DATE","DurationOfUnemployment","VALUE","DataType"]); %Extract relevant variables

df_UNEMPDUR_CAN(contains(string(df_UNEMPDUR_CAN{:,4}),'Unadjusted'),:) = []; %Remove data that is not seasonaly adjusted

df_UNEMPDUR_CAN.DurationOfUnemployment = string(df_UNEMPDUR_CAN.DurationOfUnemployment); 

df_UNEMPDUR_CAN = unstack(df_UNEMPDUR_CAN,'VALUE','DurationOfUnemployment'); %Change the data structure from long to wide

df_UNEMPDUR_CAN = df_UNEMPDUR_CAN(:,["REF_DATE","x1To4Weeks","TotalUnemployed_AllWeeks"]);%Keep only unemployment duration between 1-4 weeks and 

df_UNEMPDUR_CAN.REF_DATE = datestr(df_UNEMPDUR_CAN.REF_DATE,'yyyy-mm-dd');

df_UNEMPDUR_CAN = renamevars(df_UNEMPDUR_CAN,["REF_DATE","x1To4Weeks","TotalUnemployed_AllWeeks"],["Date","UNEMP_DUR_1_4","UNEMP_DUR_TOT"]);

% Calculate the job finding probability
df_UNEMPDUR_CAN.PROB_JOB_FINDING = NaN(height(df_UNEMPDUR_CAN),1); %Prob(A worker who starts month t unemployed finds a job within the month)
df_UNEMPDUR_CAN.OUTFLOW_RATE = NaN(height(df_UNEMPDUR_CAN),1);% The rate at which an unemployed worker finds a job

for i = 1:length(df_UNEMPDUR_CAN.Date)
   if i < length(df_UNEMPDUR_CAN.Date)
    df_UNEMPDUR_CAN.PROB_JOB_FINDING(i,1) = 1-(df_UNEMPDUR_CAN.UNEMP_DUR_TOT(i+1,1)-df_UNEMPDUR_CAN.UNEMP_DUR_1_4(i+1,1))/df_UNEMPDUR_CAN.UNEMP_DUR_TOT(i,1);
   else 
    df_UNEMPDUR_CAN.PROB_JOB_FINDING(i,1) = NaN;
   end
   df_UNEMPDUR_CAN.OUTFLOW_RATE(i,1) = (-1)*log(1-df_UNEMPDUR_CAN.PROB_JOB_FINDING(i,1));
end

%Express as percentages
df_UNEMPDUR_CAN.PROB_JOB_FINDING = df_UNEMPDUR_CAN.PROB_JOB_FINDING*100;
df_UNEMPDUR_CAN.OUTFLOW_RATE = df_UNEMPDUR_CAN.OUTFLOW_RATE*100;
%% Import Uncertainty data: EPU
df_EPU_CAN = readtable('./Data/Canada_Policy_Uncertainty_Data.xlsx','ReadRowNames', 0);

%Manipulation
df_EPU_CAN(end,:)=[]; %Remove the last row

df_EPU_CAN.Year = str2double(df_EPU_CAN.Year);

df_EPU_CAN.Date = datetime(df_EPU_CAN.Year,df_EPU_CAN.Month,1);

df_EPU_CAN = renamevars(df_EPU_CAN,"CanadaNews_BasedPolicyUncertaintyIndex","EPU_CAN");

df_EPU_CAN = df_EPU_CAN(:,["Date","EPU_CAN"]);

df_EPU_CAN.Date = datestr(df_EPU_CAN.Date,'yyyy-mm-dd');

%% Import Uncertainty data: US Macroeconomic uncertainty
df_MacroUN_US = readtable('./Data/MacroUncertaintyToCirculate.xlsx','ReadRowNames', 0);

df_MacroUN_US = renamevars(df_MacroUN_US,["Date","h_1","h_3","h_12"],["Date","US_unc_h1","US_unc_h3","US_unc_h12"]);

df_MacroUN_US.Date = datestr(df_MacroUN_US.Date, 'yyyy-mm-dd');

%% Merge tables 

% Merge the monthly data first

dfs = {df_EPU_CAN df_UNEMPDUR_CAN df_MacroUN_US};
df_merge = dfs{1};

for i = 2:numel(dfs)
    df_merge = outerjoin(df_merge, dfs{i}, 'MergeKeys',true); 
end

%%

% Fetch the variable lable
var_label = df_merge.Properties.VariableNames;

% Transform Date from Char to Datetime
df_merge.Date = datetime(df_merge.Date);

% Transform into quaterly data by first aggregated and then take the mean
df_merge = groupsummary(df_merge,"Date","quarter","mean");

% Transform the quarterly date into datetime variable
temp_time = string(df_merge.quarter_Date);

temp_time = datenum(temp_time,'QQ yyyy');

df_merge.quarter_Date = datetime(temp_time,'ConvertFrom','datenum');

df_merge.Date_str = datestr(df_merge.quarter_Date,'yyyy-mm-dd'); %Create a string for merge

%%
% Merge the quarterly data

dfs_qtr = {df_nomGDP_CAN df_realGDP_CAN df_CPI_CAN df_UNEMP_CAN df_WUI_CAN df_MacroUN_CAN};
df_merge_qtr = dfs_qtr{1};

for i = 2:numel(dfs_qtr)
    df_merge_qtr = outerjoin(df_merge_qtr, dfs_qtr{i}, 'MergeKeys',true); 
end

% Merge quaterly data with df_merge

df_merge_orig = df_merge_qtr;

df_merge_orig = outerjoin(df_merge_orig,df_merge,'MergeKeys',true);

%%
% Some Housekeeping
df_merge_orig = removevars(df_merge_orig, ["quarter_Date","GroupCount","Date_str"]);

% Export

writetable(df_merge_orig,'./Data/data_quarterly_merged.csv') 

%% Stationarity tests

var_label = df_merge_orig.Properties.VariableNames; %fetch the labels of the variables

adf_test_res = [adftest(df_merge_orig(:,2),Alpha=0.01)];

kpss_test_res = [kpsstest(df_merge_orig(:,2),Alpha=0.01)];

for j = 3:numel(var_label)
    new_adf = adftest(df_merge_orig(:,j),Alpha = 0.01);
    adf_test_res(end+1,:) = [new_adf];

    new_kpss = kpsstest(df_merge_orig(:,j),Alpha = 0.01);
    kpss_test_res(end+1,:) = [new_kpss];
end

adf_test_res.Properties.RowNames = var_label(1,2:end); %non stationary

kpss_test_res.Properties.RowNames = var_label(1,2:end); %non stationary
