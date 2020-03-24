% This script regresses TP in Maumee river with TP in Maumee bay
% Note: (1) Missing data in the river was replaced by NaN

clear all
close all
clc


direc_htlp='D:/Research/EPA_Project/Lake_Erie_HAB/Data/HTLP'; 
bay_station='4P';           % which station of maumee bay
save_direc='D:/Research/EPA_Project/Lake_Erie_HAB/matlab_codes/observed_vs_predicted_plots_03_23_2019';
datenum_wrapper=@(x)datenum(x,'mm/dd/yyyy');
frac_cal_samples=0.70;      % fraction of samples tobe used for calibration
outlier_lim=0.3;            % sample with concentration greater than this value (in bay) are considered outliers
%% read observed maumee river data
%{
fname_htlp='daily_maumeedata.txt';
filename=fullfile(direc_htlp,fname_htlp);
data=readtable(filename,'delimiter','\t');
date_river=data.(1);
datenum_river=cellfun(datenum_wrapper,date_river);
streamflow=data.Flow_CMS;
TP_conc_river=data.TP_Mg_LAsP;
TP_conc_river(TP_conc_river<0)=NaN;
TP_load_river_raw=TP_conc_river.*streamflow*(10^(-3)*24*3600);

% replace the data corresponding to missing dates with nan
max_num_datenum=[datenum_river(1):datenum_river(end)]';
TP_load_river=nan(length(max_num_datenum),1);
[~,ia,~]=intersect(max_num_datenum,datenum_river);
TP_load_river(ia)=TP_load_river_raw;
%}

%% read reconstructed maumee river data
fname='maumee_reconstructed_TP_conc.txt';
filename=fullfile(direc_htlp,fname);
data=readtable(filename,'delimiter','\t');
date_river=data.(1);
datenum_river=cellfun(datenum_wrapper,date_river);
TP_load_river=data.(4);

%% read maumee bay data
fname_lec=strcat('WQ_',bay_station);
direc_lec='D:/Research/EPA_Project/Lake_Erie_HAB/Data/lake_erie_LEC';
filename=fullfile(direc_lec,fname_lec);
data=readtable(filename,'delimiter','\t');
date_bay=data.date;
datenum_bay=cellfun(datenum_wrapper,date_bay);
TP_bay=data.tp_mg_L_;

% create predictor and response variables
empty_datenum=[];
count=0;
for date_ind=1:length(date_bay)
    
    % find datenum corresponding to bay data in river data
    datenum_tmp=datenum_bay(date_ind);
    ind=find(datenum_river==datenum_tmp);
    
    if isempty(ind)     % if no detenum couldbe found in river data
        empty_datenum=[empty_datenum;datenum_tmp];
    else
        count=count+1;
        TP_load_river_pred(:,count)=TP_load_river(ind-150:ind);       % predictor variable (each column is a predictor variable)
        TP_conc_bay(count)=TP_bay(date_ind);                        % response variable
    end
end

% remove points that correspod to NaNs in response variable
not_nan_ind=find(~isnan(TP_conc_bay));
TP_conc_bay=TP_conc_bay(not_nan_ind)';                          % response variable
TP_load_river_pred=TP_load_river_pred(:,not_nan_ind)';          % predictor variable (each row is a predictor variable)

% remove points that correspod to NaNs in predcitor variable
sum_pred=sum(TP_load_river_pred,2);
not_nan_ind=find(~isnan(sum_pred));
TP_conc_bay=TP_conc_bay(not_nan_ind);                          % response variable
TP_load_river_pred=TP_load_river_pred(not_nan_ind,:);          % predictor variable (each row is a predictor variable)

% remove outliers
TP_load_river_pred(TP_conc_bay>outlier_lim,:)=[];
TP_conc_bay(TP_conc_bay>outlier_lim)=[];


% create aggregated predictors
% pred=[TP_load_river_pred(:,1),...
%     sum(TP_load_river_pred(:,1:3),2),...
%     sum(TP_load_river_pred(:,1:6),2),...
%     sum(TP_load_river_pred(:,1:9),2),...
%     sum(TP_load_river_pred(:,1:12),2),...
%     sum(TP_load_river_pred(:,1:15),2),...
%     sum(TP_load_river_pred(:,1:18),2),...
%     sum(TP_load_river_pred(:,1:21),2),...
%     sum(TP_load_river_pred(:,1:24),2),...
%     sum(TP_load_river_pred(:,1:27),2),...
%     sum(TP_load_river_pred(:,1:30),2),...
%     sum(TP_load_river_pred(:,1:33),2),...
%     sum(TP_load_river_pred(:,1:36),2),...
%     sum(TP_load_river_pred(:,1:39),2),...
%     sum(TP_load_river_pred(:,1:42),2),...
%     sum(TP_load_river_pred(:,1:45),2)];

nsamp=length(TP_conc_bay);
cal_ind=1:round(frac_cal_samples*nsamp);      % index of samples in calibration set
val_ind=setdiff(1:nsamp,cal_ind); % rest in validation
ncal=length(cal_ind);             % number of calibration samples
nval=length(val_ind);             % number of validation samples
pred_cal=TP_load_river_pred(cal_ind,:);
resp_cal=TP_conc_bay(cal_ind);
pred_val=TP_load_river_pred(val_ind,:);
resp_val=TP_conc_bay(val_ind);
%% principle component regression
%{
[~,score,~,~,explained] = pca(TP_load_river_pred);
total_explained=cumsum(explained);
ncomp=find(total_explained>99,1,'first');

% regression analysis
mdl=fitlm(score(:,1:ncomp),TP_conc_bay);
%}



%% partial least square regression
%{
[XL,YL,XS,YS,beta] = plsregress(pred_cal,resp_cal,6);
Fitted_cal=[ones(ncal,1),pred_cal]*beta;
Fitted_val=[ones(nval,1),pred_val]*beta;
%}

%% random forest
%
NumTrees=100;
B = TreeBagger(NumTrees,pred_cal,resp_cal,...
    'Method','regression','NVarToSample',30,'MinLeaf',3); 
Fitted_cal=predict(B,pred_cal);
Fitted_val=predict(B,pred_val);
%}
%% plot
R2_cal=corr(Fitted_cal,resp_cal)^2;
R2_val=corr(Fitted_val,resp_val)^2;
R2_cal=round(100*R2_cal)/100;
R2_val=round(100*R2_val)/100;

scatter(resp_cal,Fitted_cal,'filled'); hold on
scatter(resp_val,Fitted_val,'facecolor','r','markeredgecolor','r');

axis_lim_low=min(TP_conc_bay);
axis_lim_upp=max(TP_conc_bay);
xlim([axis_lim_low axis_lim_upp])
ylim([axis_lim_low axis_lim_upp])
plot([axis_lim_low axis_lim_upp],[axis_lim_low axis_lim_upp],'color','black');
annotation('textbox', [0.15, 0.8, 0.1, 0.1], 'String',...
    {['Station: ',bay_station],['R^2_{cal} = ',num2str(R2_cal)],['R^2_{val} = ',num2str(R2_val)]},...
    'HorizontalAlignment','left','fontname','arial','fontsize',10,...
    'edgecolor','none')
xlabel('Observed TP (mg L^{-1})','fontname','arial','fontsize',10);
ylabel('Predicted TP (mg L^{-1})','fontname','arial','fontsize',10);
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',10,box);
clear box

sname=strcat('obs_vs_pred_',bay_station,'.svg');
save_filename=fullfile(save_direc,sname);
fig2svg(save_filename);