% this routine combines total CI values for composite images, obtained from 
% Sentinel data, into one text file

clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Sentinel';

product_direc = {'2016','2017','2018','2019','2020'};

datenum_wrapper = @(x)datenum(x,'yyyy-mm-dd');
CI_data=[];
datenums = [];
for dir_ind = 1:length(product_direc)
    
    filename = fullfile(direc,product_direc{dir_ind},'composite_product','total_CI.txt');
    fid = fopen(filename,'r');
    data = textscan(fid,'%s%f%f%s','delimiter','\t','headerlines',1);
    fclose(fid);
    
    CI_data = [CI_data;data{2}];
    
    dates = data{4};
    datenums = [datenums;cellfun(datenum_wrapper,dates)];
end

wfname = 'CI_total_combined_Sentinel.txt';
filename = fullfile(direc,wfname);
fid = fopen(filename,'w');
fprintf(fid,'%s\t%s\n','Begin_date(yyyy-mm-dd)','CI_total');
for wind = 1:length(CI_data)
    
    fprintf(fid,'%s\t%d\n',datestr(datenums(wind),'yyyy-mm-dd'),CI_data(wind));
    
end
fclose(fid);


