% Author: Lukas Ma (modified from codes by Andrea Bucci)
%
% Date: 2022-12-29
%
% A function for running the non-linearity test proposed by Teräsvirta & Yang (2014) 
% Based on (more like ''translated'' from) the R code written by Andrea
% Bucci. See https://github.com/andbucci/starvars/blob/master/R/VLSTARjoint.R
%


function [nonlinLMres] = nonlin_testSTVAR(data, var_par, exdata, state, state_choice, alpha)

y = data; % specify the endogenous variables
ncoly = size(y,2);

if state_choice == true
    st = state.s_multi; % specify the transition variables 
else
    st = state.s; 
end

p = var_par.p;
c_case = var_par.c_case;
alpha = alpha/100; 

VAR = estimateVAR(y, p, c_case, exdata); %VAR estimation

x = VAR.X; % Fetch the estimated X matrix
ncolx = size(x,2);
nrowx = size(x,1);


if state_choice == true
    LM3 = NaN([size(st,2) 1]);
    pvalue = NaN([size(st,2) 1]);
    st = st(VAR.p:size(st,1),:);
    for j = 1:size(st,2)
        ee = VAR.u; %reduced form residuals
        RSS0 = transpose(ee)*ee;
        ZZ = NaN([nrowx ncolx*3]);
        
        for i = 1:nrowx
            xst1 = x(i,:).*st(i,j);
            xst2 = x(i,:).*st(i,j)^2;
            xst3 = x(i,:).*st(i,j)^3;
            ZZ(i,:) = [xst1 xst2 xst3];
        end
        
        ausvar = {};
        ll = [];

        for k = 1:size(ee,2)
            temp_mod = fitlm(ZZ,ee(:,k));
            temp_res = table2array(temp_mod.Residuals(:,'Raw'));
            ausvar{k} = temp_mod;
            ll = [ll, temp_res];
        end
        
        RSS1 = transpose(ll)*ll;
        trac1 = trace(pinv(RSS0)*RSS1);
        
        LM3(j,:) = nrowx*(ncoly - trac1);
        df = 3*ncoly + (ncolx);
        conflev = 1-alpha/2;
        chi = chi2inv(conflev, df);
        pvalue(j,:) = chi2cdf(LM3(j,:), df, 'upper');
    end


else
    st = st(VAR.p:size(st,1),:);
    ee = VAR.u; %reduced form residuals
    RSS0 = transpose(ee)*ee;
    ZZ = NaN([nrowx ncolx*3]);
    for i = 1:nrowx
        xst1 = x(i,:).*st(i);
        xst2 = x(i,:).*st(i)^2;
        xst3 = x(i,:).*st(i)^3;
        ZZ(i,:) = [xst1 xst2 xst3];
    end 

    ausvar = {};
    ll = [];
    for k = 1:size(ee,2)
        temp_mod = fitlm(ZZ,ee(:,k));
        temp_res = table2array(temp_mod.Residuals(:,'Raw'));
        ausvar{k} = temp_mod;
        ll = [ll, temp_res];
    end

    RSS1 = transpose(ll)*ll;
    trac1 = trace(pinv(RSS0)*RSS1);
       
    LM3 = nrowx*(ncoly - trac1);
    df = 3*ncoly + (ncolx);
    conflev = 1-alpha/2;
    chi = chi2inv(conflev, df);
    pvalue = chi2cdf(LM3, df, 'upper');

end

% Housekeeping 
  nonlinLMres.LM3 = LM3;
  nonlinLMres.pvalue = pvalue;
  nonlinLMres.chi = chi;
  nonlinLMres.statevar = st; 
  nonlinLMres.stchoice = state_choice;
  nonlinLMres.df = df; 
  
  return
end



