% this routine plots total CI values for composite images

clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Sentinel';

product_direc = {'2016','2017','2018','2019','2020'};

xticklabels = {'20160101','20170101','20180101','20190101','20200101'};
xticklabels = cellfun(@(x)datenum(x,'yyyymmdd'),xticklabels);

split_wrapper = @(x)strsplit(x,'_');
datenum_wrapper = @(x)datenum(x,'yyyy-mm-dd');
plot_y = [];
datenums = [];
for dir_ind = 1:length(product_direc)
    
    filename = fullfile(direc,product_direc{dir_ind},'composite_product','total_CI.txt');
    fid = fopen(filename,'r');
    data = textscan(fid,'%s%f%f%s','delimiter','\t','headerlines',1);
    fclose(fid);
    
    plot_y = [plot_y;data{2}];
    
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
ylabel('Total CI','fontname','arial','fontsize',12)
box('on')
box.linewidth = 2;
set(gca,'fontname','arial','fontsize',12,'Xtick',xticklabels,'plotboxaspectratio',[2,1,1],box);
datetick('x','yyyy','keeplimits','keepticks')
clear box

sname = 'SEN_total_CI.svg';
filename = fullfile(direc,sname);
fig2svg(filename);