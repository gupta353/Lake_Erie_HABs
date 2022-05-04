%{
This script does  simpel and weighted model averaging

Author: Abhinav Gupta (Created: 4 May 2022)

%}

clear all
close all
clc

direc = 'D:/Research/EPA_project/Lake_Erie_HAB/matlab_codes/plots_04_28_2022/30_days_lead_time_prediction';

%% read LASSO, RF and ANN data
% LASSO
fname = 'LASSO_obs_pred_log_CI_cross_validation_cc10_removed_multilag_vars_included.mat';
filename = fullfile(direc, fname);
LASSO = load(filename);

% RF
fname = 'RF_obs_pred_log_CI_cross_validation_cc10_removed_multilag_vars_included.mat';
filename = fullfile(direc, fname);
RF = load(filename);

% ANN
fname = 'ANN_obs_pred_log_CI_cross_validation_cc10_removed_multilag_vars_included.mat';
filename = fullfile(direc, fname);
ANN = load(filename);

%% Averaging for bootstrap scheme
%{
CI_test_data = {};
for ii = 1:1000
    
    CI_obs = LASSO.CI_test_data{ii}(:,1);
    
    pred_LASSO = LASSO.CI_test_data{ii}(:,2);
    pred_RF = RF.CI_test_data{ii}(:,2);
    pred_ANN = ANN.CI_test_data{ii}(:,2);
    
    mse_LASSO = LASSO.cv_mse(ii);
    mse_RF = RF.cv_mse(ii);
    mse_ANN = ANN.cv_mse(ii);
    
    % simple averaging
    pred_simp_avg = (pred_LASSO + pred_RF + pred_ANN)/3;
    
    % weighted averaging
    w = [1/mse_LASSO, 1/mse_RF, 1/mse_ANN]/sum(1/mse_LASSO + 1/mse_RF + 1/mse_ANN);
    pred_wght_avg = w(1)*pred_LASSO + w(2)*pred_RF + w(3)*pred_ANN;
    
    CI_test_data{ii} = [CI_obs, pred_simp_avg, pred_wght_avg];
    
    R2_LASSO(ii) = LASSO.R21(ii);
    R2_RF(ii) = RF.R21(ii);
    R2_ANN(ii) = ANN.R21(ii);
    R2_simp(ii) = corr(CI_obs, pred_simp_avg)^2;
    R2_wght(ii) = corr(CI_obs, pred_wght_avg)^2;
    
end

% save data
fname = 'AVG_BS_data.mat';
filename = fullfile(direc, fname);
save(filename, 'R2_simp', 'R2_wght', 'CI_test_data');
%}
%% averaging for LOOCV scheme
%
Fitted_val_simp_avg = [];
Fitted_val_wght_avg = [];
for ii = 1:length(LASSO.CI)
    
    mse_LASSO = LASSO.cv_mse(ii);
    mse_RF = RF.cv_mse(ii);
    mse_ANN = ANN.cv_mse(ii);
    
    % simple averaging
    Fitted_val_simp_avg(ii,1) = (LASSO.Fitted_val(ii) + RF.Fitted_val(ii) + ANN.Fitted_val(ii))/3;
    
    
    % weighted averaging
    w = [1/mse_LASSO, 1/mse_RF, 1/mse_ANN]/sum(1/mse_LASSO + 1/mse_RF + 1/mse_ANN);
    Fitted_val_wght_avg(ii,1) = w(1)*LASSO.Fitted_val(ii) + w(2)*RF.Fitted_val(ii) + w(3)*ANN.Fitted_val(ii);
    
end

% scatter plot
CI = LASSO.CI;
Fitted_val = Fitted_val_simp_avg;
scatter(CI,Fitted_val,'filled'); hold on
ylim([2 6.5])
xlim([2 6.5])
plot([2 6.5],[2 6.5],'color','black','linewidth',2)
R2 = corr(Fitted_val,CI)^2;
R2 = (round(R2*100))/100;

xlabel('Observed log(CI)','fontname','arial','fontsize',12)
ylabel('Predicted log(CI)','fontname','arial','fontsize',12)
box('on')
box.linewidth = 2;
set(gca,'fontname','arial','fontsize',12,box)
title(['R^2 = ',num2str(R2)],'fontname','arial','fontsize',12);
clear box

% save plot
fname = 'AVG_SIMPLE_obs_pred_log_CI_cross_validation_cc10_multilag_vars_included.svg';
filename = fullfile(direc,fname);
saveas(gcf,filename,'svg')
close all
%%%
Fitted_val = Fitted_val_wght_avg;
scatter(CI,Fitted_val,'filled'); hold on
ylim([2 6.5])
xlim([2 6.5])
plot([2 6.5],[2 6.5],'color','black','linewidth',2)
R2 = corr(Fitted_val,CI)^2;
R2 = (round(R2*100))/100;

xlabel('Observed log(CI)','fontname','arial','fontsize',12)
ylabel('Predicted log(CI)','fontname','arial','fontsize',12)
box('on')
box.linewidth = 2;
set(gca,'fontname','arial','fontsize',12,box)
title(['R^2 = ',num2str(R2)],'fontname','arial','fontsize',12);
clear box

% save plot
fname = 'AVG_WEIGHT_obs_pred_log_CI_cross_validation_cc10_multilag_vars_included.svg';
filename = fullfile(direc,fname);
saveas(gcf,filename,'svg')

% save data
fname = 'AVG_SIMPLE_obs_pred_log_CI_cross_validation_cc10_multilag_vars_included.mat';
filename = fullfile(direc,fname);
save(filename, 'CI', 'Fitted_val_simp_avg', 'R2')

fname = 'AVG_WEIGHT_obs_pred_log_CI_cross_validation_cc10_multilag_vars_included.mat';
filename = fullfile(direc,fname);
save(filename, 'CI', 'Fitted_val_wght_avg', 'R2')
%}