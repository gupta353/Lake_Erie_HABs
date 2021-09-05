% script to plort probability distributions of R2 for all the different
% settings

clear all
close all
clc

direc = 'D:/Research/EPA_project/Lake_Erie_HAB/matlab_codes/plots_08_28_2021';

subdir = {'10_days_lead_time_prediction','20_days_lead_time_prediction','30_days_lead_time_prediction'};

%% read data
R2_data = [];
for sdind = 1:length(subdir)
    
    % read LASSO data
    filename = fullfile(direc,subdir{sdind},'LASSO_BS.mat');
    load(filename);
    R2_data = [R2_data,R21'];
    
    % read RF data
    filename = fullfile(direc,subdir{sdind},'RF_BS.mat');
    load(filename);
    R2_data = [R2_data,R21'];
   
end

% boxplot
xvar = ['LASSO';'RF   ';'LASSO';'RF   ';'LASSO';'RF   '];
boxplot(R2_data,'Labels',xvar);
box('on')
box.linewidth = 2;
set(gca,'fontname','arial','fontsize',10,box);
xlabel('Regression method','fontname','arial','fontsize',10);
ylabel('R^2','fontname','arial','fontsize',10);

% save plot
sname = 'R2_BS_combined.svg';
filename = fullfile(direc,sname);
fig2svg(filename);