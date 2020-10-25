% statitical model to estimate 10-day timescale CI

clear all
close all
clc

% read data
fname = 'model_data.txt';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB','matlab_codes',fname);
fid = fopen(filename,'r');
formatspec = ['%s',repmat('%f',1,18)];
data = textscan(fid,formatspec,'delimiter','\t','headerlines',1);
fclose(fid);

dates = data{1};
CI = data{2};
preds = cat(2,data{3:end});
preds(:,1:6) = [];
wrapper = @(x)str2num(datestr(datenum(x,'dd-mmm-yyyy'),'mm'));
month_num = cellfun(wrapper,dates);
% preds = [preds,month_num];
pred_name = {'CI','Minimum wind speed','Maximum air temperature','Average TP','Average TKN','Ratio of TP to TKN','Average streamflow','Secchi depth','TP load from Mar to Jun','TKN load from Jan to Jun'};


% read cloud cover data
fname = 'cloud_cover_combined_MERIS.txt';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB','Data/remote_sensing_data',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid)
cc_dates = data{1};
cc = data{2};

% remove outliers and missing values
zero_ind = find(cc>0.10 | isnan(preds(:,2)));
CI(zero_ind) = [];
preds(zero_ind,:) = [];
dates(zero_ind) = [];

% normalize the predictor variables
preds_m = mean(preds);
preds_std = std(preds);
preds = bsxfun(@minus,preds,preds_m);
preds = preds./repmat(preds_std,size(preds,1),1);
preds(:,end) = [];

% calibration and validation data
%{
rng(1);
count = 0;
while count<1000
    count = count+1;
    cal_ind = randsample(length(CI),60);
    val_ind = setdiff(1:length(CI),cal_ind);
    CI_cal = log(CI(cal_ind)); preds_cal = preds(cal_ind,:);
    CI_val = log(CI(val_ind)); preds_val = preds(val_ind,:);
    
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
    [B,Fitinfo] = lasso(preds_cal,CI_cal,'alpha',0.0001,'CV',10);
    ind = find(Fitinfo.MSE == min(Fitinfo.MSE));
    beta(:,count) = [Fitinfo.Intercept(ind);B(:,ind)];
    Fitted_val = beta(1,count) + preds_val*beta(2:end,count);
    
%     scatter(CI_val,Fitted_val,'filled'); hold on
%     ylim([0 10])
%     xlim([0 10])
%     plot([0 10],[0 10],'color','black','linewidth',2)
%     
%     R21(count)=corr(CI_val,Fitted_val)^2;
%     lambda(count) = Fitinfo.LambdaMinMSE;
%     R21_tmp = (round(R21(count)*100))/100;
%     title(['R^2 = ',num2str(R21_tmp)]);
%     
%     xlabel('Observed log(CI)','fontname','arial','fontsize',12)
%     ylabel('Predicted log(CI)','fontname','arial','fontsize',12)
%     box('on')
%     box.linewidth = 2;
%     set(gca,'fontname','arial','fontsize',12,box)
%     pause;
%     clear box
    
    %     fname = strcat('obs_pred_log_CI_',num2str(count),'.jpg');
    %     filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots/ridge_regression_plots_2',fname);
    %     print(filename,'-r300','-djpeg')
%     close all
    
end
%}
% hist(R21)
% xlabel('Coefficient of determination (R^2)','fontname','arial','fontsize',12);
% ylabel('Number of samples in the bin','fontname','arial','fontsize',12)
% fname = strcat('R2_histogram.svg');
% filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots/ridge_regression_plots_2',fname);
% fig2svg(filename);
% %
% hist(lambda)
% xlabel('Regularization parameter (\lambda)','fontname','arial','fontsize',12);
% ylabel('Number of samples in the bin','fontname','arial','fontsize',12)
% fname = strcat('lambda_histogram.svg');
% filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots/ridge_regression_plots_2',fname);
% fig2svg(filename);
%{
for pind = 1:length(pred_name)
    
    hist(beta(1+pind,:))
    xlabel(pred_name{pind},'fontname','arial','fontsize',12);
end

boxplot(beta(2:end,:)',pred_name);
ylabel('Regression coefficients','fontname','arial','fontsize',12)
fname = 'coefficient_distribtuion.svg';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/plots/ridge_regression_plots_2',fname);
fig2svg(filename);
%}
%}
%% cross-validation (no validation) lasso
%
CI = log(CI);
for CI_ind =  1:length(CI)
    
    val_ind = CI_ind;
    cal_ind = setdiff(1:length(CI),val_ind);
    CI_cal = CI(cal_ind); preds_cal = preds(cal_ind,:);
    CI_val = CI(val_ind); preds_val = preds(val_ind,:);
    
    [B,Fitinfo] = lasso(preds_cal,CI_cal,'alpha',0.0001,'CV',10);
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