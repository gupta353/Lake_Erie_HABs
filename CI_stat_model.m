% statitical model to estimate 10-day timescale CI

clear all
close all
clc

% read data
fname = 'model_data.txt';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB','matlab_codes',fname);
fid = fopen(filename,'r');
formatspec = ['%s',repmat('%f',1,63)];
data = textscan(fid,formatspec,'delimiter','\t','headerlines',1);
fclose(fid);

dates = data{1};
CI = data{2};
preds = cat(2,data{3:end});
preds(:,1:24) = [];
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

% covert CIs into logarithmic scale
CI = log(CI);
preds(:,2) = log(preds(:,2));

% remove outliers and missing values
zero_ind = find(cc>0.10 | isnan(preds(:,2)));
CI(zero_ind) = [];
preds(zero_ind,:) = [];
dates(zero_ind) = [];

% normalize the predictor variables
preds_m = nanmean(preds);
preds_std = nanstd(preds);
preds = bsxfun(@minus,preds,preds_m);
preds = preds./repmat(preds_std,size(preds,1),1);
% preds(:,end) = [];

% calibration and validation data
%{
rng(1);
count = 0;
for count = 1:1000
    
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
    %{
% cross-validation
NumTrees=50:80;
NVarToSample=2:15;          % number of predictors that random forest considers at each node
MinLeaf=1:5;
count = 0;

for par1_ind = 1:length(NumTrees)
    for par2_ind = 1:length(NVarToSample)
        for par3_ind = 1:length(MinLeaf)
            count = count+1;
            
            %   predfun = @(Xcal,ycal,Xval)RFreg(Xcal,ycal,Xval,NumTrees(par1_ind),NVarToSample(par2_ind),MinLeaf(par3_ind));
            %   err = crossval('mse',preds_cal,CI_cal,'predfun',predfun);
            %   err_mat(count,:) = [NumTrees(par1_ind),NVarToSample(par2_ind),MinLeaf(par3_ind),err];
            B = TreeBagger(NumTrees(par1_ind),preds_cal,CI_cal,...
                'Method','regression','NVarToSample',NVarToSample(par2_ind),'MinLeaf',MinLeaf(par3_ind));
            Fitted_val=predict(B,preds_val);
            scatter(CI_val,Fitted_val); hold on
            ylim([0 10])
            xlim([0 10])
            plot([0 10],[0 10],'color','black')
            R2(count,:)=[NumTrees(par1_ind),NVarToSample(par2_ind),MinLeaf(par3_ind),corr(CI_val,Fitted_val)^2];
%             title(['R^2 = ',num2str(R2(count,4))])
%             pause;
%             close all
        end
    end
end

ind = find(R2(:,4) == max(R2(:,4)));
B = TreeBagger(R2(ind,1),preds_cal,CI_cal,...
    'Method','regression','NVarToSample',R2(ind,2),'MinLeaf',R2(ind,3));
Fitted_val=predict(B,preds_val);
scatter(CI_val,Fitted_val)
ylim([0 700])
xlim([0 700])
R21=corr(CI_val,Fitted_val)^2;
title(['R^2 = ',num2str(R21)])
    %}
    
    % linear regression
    %{
mdl=fitlm(preds_cal,CI_cal);
Fitted_val = predict(mdl,preds_val);

scatter(CI_val,Fitted_val); hold on
ylim([0 10])
xlim([0 10])
plot([0 10],[0 10],'color','black')

R21=corr(CI_val,Fitted_val)^2;
    %}
    
    % lasso regression
    %
    [B,Fitinfo] = lasso(preds_cal,CI_cal,'alpha',0.999,'CV',5);
    ind = find(Fitinfo.MSE == min(Fitinfo.MSE));
    beta(:,count) = [Fitinfo.Intercept(ind);B(:,ind)];
    Fitted_val =  [ones(size(preds_val,1),1) preds_val]*beta(:,count);
    
    scatter(CI_val,Fitted_val,'filled'); hold on
    ylim([2 8])
    xlim([2 8])
    plot([0 10],[0 10],'color','black','linewidth',2)
    
    R21(count)=corr(CI_val,Fitted_val)^2;
    lambda(count) = Fitinfo.LambdaMinMSE;
    R21_tmp = (round(R21(count)*100))/100;
    title(['R^2 = ',num2str(R21_tmp)]);
    
    xlabel('Observed log(CI)','fontname','arial','fontsize',12)
    ylabel('Predicted log(CI)','fontname','arial','fontsize',12)
    box('on')
    box.linewidth = 2;
    set(gca,'fontname','arial','fontsize',12,box)
%     pause(1);
    clear box
    
    fname = strcat('obs_pred_log_CI_',num2str(count),'.jpg');
    filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots/lasso_regression_plots_1',fname);
    print(filename,'-r300','-djpeg')
    close all
    
end
%
hist(R21)
xlabel('Coefficient of determination (R^2)','fontname','arial','fontsize',12);
ylabel('Number of samples in the bin','fontname','arial','fontsize',12)
fname = strcat('R2_histogram.svg');
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots/lasso_regression_plots_1',fname);
fig2svg(filename);
%
hist(lambda)
xlabel('Regularization parameter (\lambda)','fontname','arial','fontsize',12);
ylabel('Number of samples in the bin','fontname','arial','fontsize',12)
fname = strcat('lambda_histogram.svg');
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots/lasso_regression_plots_1',fname);
fig2svg(filename);
%
for pind = 1:length(pred_name)
    
    hist(beta(1+pind,:))
    xlabel(num2str(pind),'fontname','arial','fontsize',12);
end

boxplot(beta(2:end,:)',1:38);
ylabel('Regression coefficients','fontname','arial','fontsize',12)
set(gca,'fontname','arial','fontsize',10,'plotboxaspectratio',[2 1 1])
fname = 'coefficient_distribtuion.svg';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots/lasso_regression_plots_1',fname);
fig2svg(filename);
%}
%}
%% cross-validation (no validation) lasso
%
for CI_ind =  1:length(CI)
    
    val_ind = CI_ind;
    cal_ind = setdiff(1:length(CI),val_ind);
    CI_cal = CI(cal_ind); preds_cal = preds(cal_ind,:);
    CI_val = CI(val_ind); preds_val = preds(val_ind,:);
    
    [B,Fitinfo] = lasso(preds_cal,CI_cal,'alpha',0.999,'CV',5);
    ind = find(Fitinfo.MSE == min(Fitinfo.MSE));
    beta = [Fitinfo.Intercept(ind);B(:,ind)];
    Fitted_val(CI_ind,1) = beta(1) + preds_val*beta(2:end);
    
end

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

fname = 'obs_pred_log_CI_cross_validation_cc10_removed.svg';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots',fname);
fig2svg(filename)
%}
% Uncertainty analysis of residuals
%
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

fname = 'lasso_uncertainty_analysis_cross_validation_cc10_removed.svg';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots',fname);
fig2svg(filename)

% QQ plot
figure;
qqplot(res);
box('on')
box.linewidth = 2;
set(gca,'fontname','arial','fontsize',12,box)
clear box

fname = 'qqplot_cross_validation_cc10_removed.svg';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots',fname);
fig2svg(filename)
%}


