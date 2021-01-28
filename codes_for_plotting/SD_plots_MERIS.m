% this routine plots total CI values for composite images

clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data';

product_direc = {'gupta353_MERIS_full_resolution_L2_2002_001_2020-05-21T14-50-31',...
    'gupta353_MERIS_full_resolution_L2_2003_001_2020-05-21T14-58-13',...
    'gupta353_MERIS_full_resolution_L2_2004_001_2020-05-21T15-42-46',...
    'gupta353_MERIS_full_resolution_L2_2005_001_2020-05-19T20-05-41',...
    'gupta353_MERIS_full_resolution_L2_2006_001_2020-05-19T22-11-38',...
    'gupta353_MERIS_full_resolution_L2_2007_001_2020-05-19T23-52-06',...
    'gupta353_MERIS_full_resolution_L2_2008_001_2020-05-20T01-12-24'...
    'gupta353_MERIS_full_resolution_L2_2009_001_2020-05-20T23-31-17'...
    'gupta353_MERIS_full_resolution_L2_2010_001_2020-05-21T00-09-03'...
    'gupta353_MERIS_full_resolution_L2_2011_001_2020-05-21T00-51-46'};

xticklabels = {'20020101','20030101','20040101','20050101','20060101','20070101','20080101','20090101','20100101','20110101'};
xticklabels = cellfun(@(x)datenum(x,'yyyymmdd'),xticklabels);

split_wrapper = @(x)strsplit(x,'_');
datenum_wrapper = @(x)datenum(x,'yyyy-mm-dd');
plot_y = [];
datenums = [];
for dir_ind = 1:length(product_direc)
    
    filename = fullfile(direc,product_direc{dir_ind},'composite_sd_product_1','total_secchi_depth_over_pos_CI.txt');
    fid = fopen(filename,'r');
    data = textscan(fid,'%s%f%f%s%f','delimiter','\t','headerlines',1);
    fclose(fid);
    
    plot_y = [plot_y;data{5}];
    
    dates = data{4};
    datenums = [datenums;cellfun(datenum_wrapper,dates)];
end

scatter(datenums,plot_y,5,'filled')
hold on
for date_ind = 1:length(xticklabels)
    plot([xticklabels(date_ind) xticklabels(date_ind)],[0 max(plot_y)],'--','color','black')
end
ylim([0 max(plot_y)])
xlabel('Year (first day of the year is marked)','fontname','arial','fontsize',12)
ylabel('Total Secchi disk depth','fontname','arial','fontsize',12)
box('on')
box.linewidth = 2;
set(gca,'fontname','arial','fontsize',12,'Xtick',xticklabels,'plotboxaspectratio',[2,1,1],box);
datetick('x','yyyy','keeplimits','keepticks')
clear box

sname = 'MERIS_total_SD.svg';
filename = fullfile(direc,sname);
fig2svg(filename);