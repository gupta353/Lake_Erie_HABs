% this routine combines average secchi depth values for composite images, obtained form
% Sentinel data, into one text file

clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Sentinel';

product_direc = {'2016','2017','2018','2019','2020'};

datenum_wrapper = @(x)datenum(x,'yyyy-mm-dd');
avg_SD_data=[];
datenums = [];
for dir_ind = 1:length(product_direc)
    
    filename = fullfile(direc,product_direc{dir_ind},'composite_sd_product','total_secchi_depth_over_pos_CI.txt');
    fid = fopen(filename,'r');
    data = textscan(fid,'%s%f%f%s%f','delimiter','\t','headerlines',1);
    fclose(fid);
    
    avg_SD_data = [avg_SD_data;data{5}];
    
    dates = data{4};
    datenums = [datenums;cellfun(datenum_wrapper,dates)];
    
end

wfname = 'secchi_depth_over_pos_CI_total_combined_SEN.txt';
filename = fullfile(direc,wfname);
fid = fopen(filename,'w');
fprintf(fid,'%s\t%s\n','Begin_date(yyyy-mm-dd)','average_secchi_depth(m)');
for wind = 1:length(avg_SD_data)
    
    fprintf(fid,'%s\t%d\n',datestr(datenums(wind),'yyyy-mm-dd'),avg_SD_data(wind));
    
end
fclose(fid);


