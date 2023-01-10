% Codes to produce the results presented in 
% Term paper EC7415: Uncertainty shocks and unemployment dynamics in Canada
%
% Current version: 2023-01-10
% -----------------------------------------------------------------------
% Uncomment to run a particular script. The default setting will only run
% the data preparation, transformation, and the baseline analysis 1
% -----------------------------------------------------------------------
% The cleaned data is already created in the Data folder. One can skip the
% data preparation scripts if he or she so wish. 
%

%--------------------------------------------------------------------------
%Data preparation
%--------------------------------------------------------------------------

data_prep_raw_quarter % Merge raw data from different sources
data_transform_quarter % Tranform the data

%--------------------------------------------------------------------------
%Statistical tests
%--------------------------------------------------------------------------

% Non-linearity tests

%multivar_nonlinear_test_1 %Perform nonlinearity test for baseline model 1
%multivar_nonlinear_test_2 %Perform nonlinearity test for baseline model 2

% Stationarity tests

%stationarity_tests 

%--------------------------------------------------------------------------
%Plots
%--------------------------------------------------------------------------

%figure_1_recession %Plot Figure 1 of the paper

%--------------------------------------------------------------------------
% Data analysis (baseline models)
% --
% The default setting is a 1 StD shock to US uncertainty, to evaluate a
% shock to Canadian uncertainty, use "shockpos = 5;" in the the
% Parameterization section of the script.
%--------------------------------------------------------------------------

stvar_canada_baseline_1 %Estimate a VAR and STVAR according to Eq.(5). Default setting: US uncertainty shock;

%stvar_canada_baseline_2 %Estimate a VAR and STVAR according to Eq.(6). Default setting: US uncertainty shock

%--------------------------------------------------------------------------
%Robustness checks
%--------------------------------------------------------------------------

%---------- Alternative ordering ----------------------%

%stvar_canada_robust_1 %Default setting: Canadian uncertainty shock
%stvar_canada_robust_2 %Default setting: Canadian uncertainty shock

%---------- Alternative uncertainty measure -----------%

%stvar_canada_robust_3 % Robustness check using one-year-ahead uncertainty. Default setting: CAN uncertainty shock

%---------- Alternative transition variable ------------%

%stvar_canada_robust_4 % Robustness check using a MA(4) of real GDP growth. Default setting: CAN uncertainty shock

