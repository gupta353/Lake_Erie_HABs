% statitical model to estimate 10-day timescale CI

clear all
close all
clc

% read data
fname = 'model_data_10_corr_CI_lag_included_and_TP_TKN_previous_120_and_time_step_included.txt';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB','matlab_codes',fname);
fid = fopen(filename,'r');
formatspec = ['%s',repmat('%f',1,82)];
data = textscan(fid,formatspec,'delimiter','\t','headerlines',1);
fclose(fid);

dates = data{1};
CI = data{2};
preds = cat(2,data{3:end});
preds(:,1:24) = [];
preds(:,27:34) = []; % remove spring TP and TKN loads
preds(:,31:40) = []; % remove correlation-lag variables
% preds(:,end)=[];     % remove time-step of the 10-day time-period window
wrapper = @(x)str2num(datestr(datenum(x,'dd-mmm-yyyy'),'mm'));
month_num = cellfun(wrapper,dates);
% preds = [preds,month_num];
pred_name = {'CI','Minimum wind speed','Average wind speed','Maximum wind speed','Minimum air temperature','Average air temperature','Maximum air temperature',...
    'Average TP Maumee','Average TP Raisin','Average TP Sandusky','Average TP Cuyahoga',...
    'Average TKN Maumee','Average TKN Raisin','Average TKN Sandusky','Average TKN Cuyahoga',...
    'Average streamflow Maumee','Average streamflow Raisin','Average streamflow Sandusky','Average streamflow Cuyahoga',...
    'Ratio of TP to TKN Maumee','Ratio of TP to TKN Raisin','Ratio of TP to TKN Sandusky','Ratio of TP to TKN Cuyahoga',...
    'Average solar radiation','Average water level','Secchi depth',...
    'Spring TP Maumee','Spring TP Raisin','Spring TP Sandusky','Spring TP Cuyahoga',...
    'Jan to Jun TKN Maumee','Jan to Jun TKN Raisin','Jan to Jun TKN Sandusky','Jan to Jun TKN Cuyahoga',...
    'Legacy TP Maumee','Legacy TP Raisin','Legacy TP Sandusky','Legacy TP Cuyahoga'};

% read cloud cover data
fname = 'cloud_cover_MERIS_SEN.txt';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB','Data/remote_sensing_data',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid)
cc_dates = data{1};
cc = data{2};
% remove 2020 data
cc_datenums = cellfun(@(x)datenum(x,'yyyy-mm-dd'),cc_dates);
ind = find(cc_datenums>=datenum('01/01/2020','mm/dd/yyyy'));
cc_dates(ind)  = [];
cc(ind) = [];

% remove outliers, zero CI values and missing values
zero_ind = find(cc>0.10 | isnan(preds(:,2)) | preds(:,1)==0 | isnan(preds(:,26)) | isnan(preds(:,end)) | preds(:,end)==0);
CI(zero_ind) = [];
preds(zero_ind,:) = [];
dates(zero_ind) = [];

% covert CIs into logarithmic scale
CI = log(CI);
% preds(:,1) = log(preds(:,1));
% preds(:,2:4) = preds(:,2:4).^2;
% preds = [preds(:,1) preds(:,end)];

% normalize the predictor variables
preds_m = nanmean(preds);
preds_std = nanstd(preds);
preds = bsxfun(@minus,preds,preds_m);
preds = preds./repmat(preds_std,size(preds,1),1);
% preds(:,end) = [];

% calibration and validation data
%
rng(1);

for count = 1:1000
    
    disp(count)
    cal_ind = randsample(length(CI),107);
    val_ind = setdiff(1:length(CI),cal_ind);
    CI_cal = CI(cal_ind); preds_cal = preds(cal_ind,:);
    CI_val = CI(val_ind); preds_val = preds(val_ind,:);
    
    % stratified sampling
    % data = [CI,preds];
    % data = sortrows(data);
    % val_ind = 2:3:size(data,1);
    % cal_ind = setdiff(1:length(CI),val_ind);
    % CI_cal = log(CI(cal_ind)); preds_cal = preds(cal_ind,:);
    % CI_val = log(CI(val_ind)); preds_val = preds(val_ind,:);
    %% random forest
%
    % cross-validation
    NumTrees=25:25:100;
    NVarToSample=4:4:16;          % number of predictors that random forest considers at each node
    MinLeaf=2:2:6;
    
    rcount = 0;
    for par1_ind = 1:length(NumTrees)
        for par2_ind = 1:length(NVarToSample)
            for par3_ind = 1:length(MinLeaf)
                rcount = rcount+1;
                
                %   predfun = @(Xcal,ycal,Xval)RFreg(Xcal,ycal,Xval,NumTrees(par1_ind),NVarToSample(par2_ind),MinLeaf(par3_ind));
                %   err = crossval('mse',preds_cal,CI_cal,'predfun',predfun);
                %   err_mat(count,:) = [NumTrees(par1_ind),NVarToSample(par2_ind),MinLeaf(par3_ind),err];
%                 B = TreeBagger(NumTrees(par1_ind),preds_cal,CI_cal,...
%                     'Method','regression','NVarToSample',NVarToSample(par2_ind),'MinLeaf',MinLeaf(par3_ind));
%                 Fitted_val=predict(B,preds_val);
%                 scatter(CI_val,Fitted_val); hold on
%                 ylim([0 10])
%                 xlim([0 10])
%                 plot([0 10],[0 10],'color','black')
%                 R2(rcount,:)=[NumTrees(par1_ind),NVarToSample(par2_ind),MinLeaf(par3_ind),corr(CI_val,Fitted_val)^2];
                %             title(['R^2 = ',num2str(R2(count,4))])
                %             pause;
                %             close all
                predfun = @(Xtrain,ytrain,Xtest)RFag(Xtrain,ytrain,Xtest,NumTrees(par1_ind),NVarToSample(par2_ind),MinLeaf(par3_ind));
                val = crossval('mse',preds_cal,CI_cal,'predfun',predfun);
                mse(rcount,:) = [NumTrees(par1_ind),NVarToSample(par2_ind),MinLeaf(par3_ind),val];
            end
        end
    end
    
    ind = find(mse(:,4) == min(mse(:,4)));
    B = TreeBagger(mse(ind,1),preds_cal,CI_cal,...
        'Method','regression','NVarToSample',mse(ind,2),'MinLeaf',mse(ind,3),'oobvarimp','on');
    Fitted_val=predict(B,preds_val);
    importance(count,:) = B.OOBPermutedVarDeltaError;
    scatter(CI_val,Fitted_val,'filled'); hold on
    ylim([2 8])
    xlim([2 8])
    plot([0 10],[0 10],'color','black','linewidth',2); hold off
    
    R21(count) = corr(CI_val,Fitted_val)^2;
    R21(count) = round(R21(count)*100)/100;
    title(['R^2 = ',num2str(R21(count))])
    
    xlabel('Observed log(CI)','fontname','arial','fontsize',12)
    ylabel('Predicted log(CI)','fontname','arial','fontsize',12)
    box('on');
    box.linewidth = 2;
    set(gca,'fontname','arial','fontsize',12,box)
    
    % save figure
    fname = strcat('obs_pred_log_CI_',num2str(count),'.jpg');
    filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots_08_28_2021/RF_regression_plots_BS',fname);
    print(filename,'-r300','-djpeg')
    clear box
    
%}
    
    %% linear regression
%{
    mdl=fitlm(preds_cal,CI_cal);
    Fitted_val = predict(mdl,preds_val);
    
    scatter(CI_val,Fitted_val); hold on
    ylim([0 10])
    xlim([0 10])
    plot([0 10],[0 10],'color','black')
    
    R21=corr(CI_val,Fitted_val)^2;
%}
    
    %% lasso regression
%{
    [B,Fitinfo] = lasso(preds_cal,CI_cal,'alpha',0.999,'CV',5);
    ind = find(Fitinfo.MSE == min(Fitinfo.MSE));
    beta(:,count) = [Fitinfo.Intercept(ind);B(:,ind)];
    Fitted_val =  [ones(size(preds_val,1),1) preds_val]*beta(:,count);
    
    scatter(CI_val,Fitted_val,'filled'); hold on
    ylim([2 8])
    xlim([2 8])
    plot([0 10],[0 10],'color','black','linewidth',2)
    
    R21(count) = corr(CI_val,Fitted_val)^2;
    rmse(count) = (mean((CI_val-Fitted_val).^2))^0.5;
    lambda(count) = Fitinfo.LambdaMinMSE;
    
    R21_tmp = (round(R21(count)*100))/100;
    title(['R^2 = ',num2str(R21_tmp)]);

    xlabel('Observed log(CI)','fontname','arial','fontsize',12)
    ylabel('Predicted log(CI)','fontname','arial','fontsize',12)
    box('on');
    box.linewidth = 2;
    set(gca,'fontname','arial','fontsize',12,box)
%     pause(1);
    clear box

    fname = strcat('obs_pred_log_CI_',num2str(count),'.jpg');
    filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots_08_28_2021/lasso_regression_plots_BS',fname);
    print(filename,'-r300','-djpeg')
    close all
%}
    
    %% Gaussian processes
    % divide the entire data into k uniform folds
%{
    k = 1;                                                % number of folds for cross-validation
    c = cvpartition(CI_cal,'holdout',0.5);
    % minimize objective function
    loss=@(theta)GPRobj(theta,preds_cal,CI_cal,k,c);      % loss function
    parent = ones(1,40);
    lb = 0.001*ones(1,40);                 % lower bound
    ub = 1000*ones(1,40);                                  % upper bound
    options = optimset('TolFun',10^-8,'MaxFunEvals',1000);
    [theta_opt,fval] = simulannealbnd(loss,parent,lb,ub,options);
    
    train_idx = training(c,1);
    test_idx = test(c,1);
    Fitted_val = GPPred(theta_opt,preds_cal(train_idx,:),CI_cal(train_idx,:),preds_val)';
    
    scatter(CI_val,Fitted_val,'filled'); hold on
    ylim([2 8])
    xlim([2 8])
    plot([0 10],[0 10],'color','black','linewidth',2); hold off
    
    R21(count)=corr(CI_val,Fitted_val)^2;
    R21_tmp = (round(R21(count)*100))/100;
    title(['R^2 = ',num2str(R21_tmp)]);
    
    xlabel('Observed log(CI)','fontname','arial','fontsize',12)
    ylabel('Predicted log(CI)','fontname','arial','fontsize',12)
    box('on')
    box.linewidth = 2;
    set(gca,'fontname','arial','fontsize',12,box)
    pause;
    clear box
    
end
%}
end
%
hist(R21)
xlabel('Coefficient of determination (R^2)','fontname','arial','fontsize',12);
ylabel('Number of samples in the bin','fontname','arial','fontsize',12)
fname = strcat('R2_histogram.svg');
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots_08_28_2021/RF_regression_plots_BS',fname);
saveas(gcf,filename,'svg');
%
% hist(lambda)
% xlabel('Regularization parameter (\lambda)','fontname','arial','fontsize',12);
% ylabel('Number of samples in the bin','fontname','arial','fontsize',12)
% fname = strcat('lambda_histogram.svg');
% filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots_08_28_2021/lasso_regression_plots_BS',fname);
% saveas(gcf,filename,'svg');
% %
% for pind = 1:length(pred_name)
%
%     hist(beta(1+pind,:))
%     xlabel(num2str(pind),'fontname','arial','fontsize',12);
% end

% beta_plot = beta(2:end,:);
% boxplot(beta_plot',1:39);
% ylabel('Regression coefficients','fontname','arial','fontsize',12)
% set(gca,'fontname','arial','fontsize',10,'plotboxaspectratio',[2 1 1])
% fname = 'coefficient_distribtuion.svg';
% filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots_08_28_2021/lasso_regression_plots_BS',fname);
% saveas(gcf,filename,'svg');

boxplot(importance,1:39);
ylabel('Predictor importance','fontname','arial','fontsize',12)
set(gca,'fontname','arial','fontsize',10,'plotboxaspectratio',[2 1 1])
fname = 'predictor_importance.svg';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots_08_28_2021/RF_regression_plots_BS',fname);
saveas(gcf,filename,'svg');
%}
%}
%% cross-validation (no validation)
%{
for CI_ind =  1:length(CI)
    CI_ind
    val_ind = CI_ind;
    cal_ind = setdiff(1:length(CI),val_ind);
    CI_cal = CI(cal_ind); preds_cal = preds(cal_ind,:);
    CI_val = CI(val_ind); preds_val = preds(val_ind,:);
    
    %% LASSO
%{
    [B,Fitinfo] = lasso(preds_cal,CI_cal,'alpha',0.999,'CV',5);
    ind = find(Fitinfo.MSE == min(Fitinfo.MSE));
    beta = [Fitinfo.Intercept(ind);B(:,ind)];
    Fitted_val(CI_ind,1) = beta(1) + preds_val*beta(2:end);
%}
    
    %% Random Forest
%
     NumTrees=25:25:100;
    NVarToSample=4:4:16;          % number of predictors that random forest considers at each node
    MinLeaf=2:2:6;
    
    rcount = 0;
    for par1_ind = 1:length(NumTrees)
        for par2_ind = 1:length(NVarToSample)
            for par3_ind = 1:length(MinLeaf)
                rcount = rcount+1;
                
                predfun = @(Xtrain,ytrain,Xtest)RFag(Xtrain,ytrain,Xtest,NumTrees(par1_ind),NVarToSample(par2_ind),MinLeaf(par3_ind));
                val = crossval('mse',preds_cal,CI_cal,'predfun',predfun);
                mse(rcount,:) = [NumTrees(par1_ind),NVarToSample(par2_ind),MinLeaf(par3_ind),val];
            
            end
        end
    end
    
    ind = find(mse(:,4) == min(mse(:,4)));
    B = TreeBagger(mse(ind,1),preds_cal,CI_cal,...
        'Method','regression','NVarToSample',mse(ind,2),'MinLeaf',mse(ind,3),'oobvarimp','on');
    Fitted_val(CI_ind,1)=predict(B,preds_val);

%
end
%}

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
fname = 'RF_obs_pred_log_CI_cross_validation_cc10_removed_time_steps_predictor_added.svg';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots_08_28_2021',fname);
saveas(gcf,filename,'svg')

% save data
fname = 'RF_obs_pred_log_CI_cross_validation_cc10_removed_time_steps_predictor_added.mat';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots_08_28_2021',fname);
save(filename);
%}
% Uncertainty analysis of residuals
%{
res = CI - Fitted_val;
std_dev = std(res);
upper_bound = Fitted_val + 2*std_dev;
lower_bound = Fitted_val - 2*std_dev;

% uncertainty plot
figure;
plot(Fitted_val,'linewidth',2); hold on;
scatter(1:length(CI),CI,'filled')
plot(upper_bound,'r','linewidth',2)
plot(lower_bound,'r','linewidth',2)

xlabel('Sample number','fontname','arial','fontsize',12)
ylabel('log(CI)','fontname','arial','fontsize',12)
box('on')
box.linewidth = 2;
set(gca,'fontname','arial','fontsize',12,'plotboxaspectratio',[2 1 1],'xlim',[0 180],box)
clear box
legend({'Predicted','Observed','95% prediction interval'},'fontname','arial','fontsize',12,'location','northwest')
legend('boxoff');

fname = 'RF_uncertainty_analysis_cross_validation_cc10_removed.svg';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots',fname);
fig2svg(filename)

% QQ plot
figure;
qqplot(res);
box('on')
box.linewidth = 2;
set(gca,'fontname','arial','fontsize',12,box)
clear box

fname = 'RF_qqplot_cross_validation_cc10_removed.svg';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots',fname);
fig2svg(filename)
%}

%% regression using important predictor variables
%{
% read beta values obtained by LASSO
fname = 'LASSO_BS_corr_lag_CI_beta_values.mat';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots','10_days_ahead_of_time',fname);
beta_vals = load(filename);
beta_vals = beta_vals.beta;
% remove bias term
beta_vals(1,:) = [];

% compute median of beta values
beta_M = mean(beta_vals,2);
beta_M = beta_M/max(beta_M);
beta_M = abs(beta_M);
beta_M = [beta_M,[1:48]'];
beta_M = sortrows(beta_M);
beta_M = flipud(beta_M);

%{
for var_ind = 1:size(beta_M,1)
    preds_tmp = preds(:,beta_M(1:var_ind,2));
    
    for CI_ind =  1:length(CI)
        
        val_ind = CI_ind;
        cal_ind = setdiff(1:length(CI),val_ind);
        CI_cal = CI(cal_ind); preds_cal = preds_tmp(cal_ind,:);
        CI_val = CI(val_ind); preds_val = preds_tmp(val_ind,:);
        
        %% LASSO
        %
        [B,Fitinfo] = lasso(preds_cal,CI_cal,'alpha',0.999,'CV',5);
        ind = find(Fitinfo.MSE == min(Fitinfo.MSE));
        beta = [Fitinfo.Intercept(ind);B(:,ind)];
        Fitted_val(CI_ind,1) = beta(1) + preds_val*beta(2:end);
        %}
    end
    mse(var_ind) = mean((CI-Fitted_val).^2);
    R21(var_ind) = corr(CI,Fitted_val)^2;
end

% plot mean-square-error 
plot(mse,'-o','linewidth',1,'color','black','Markeredgecolor','b','Markerfacecolor','b');
xlabel('Predictor variable','fontname','arial','fontsize',12)
ylabel('Mean-squared-error','fontname','arial','fontsize',12)
box('on')
box.linewidth = 2;
set(gca,'fontname','arial','fontsize',12,'plotboxaspectratio',[2 1 1],'ylim',[min(mse)-0.001 max(mse)+0.001],box)
clear box

% save plot
sname = 'mse_predictor_importance.svg';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots','10_days_ahead_of_time',sname);
fig2svg(filename);

% save data
sname = 'mse_predictor_importance.mat';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots','10_days_ahead_of_time',sname);
save(filename);

%}


