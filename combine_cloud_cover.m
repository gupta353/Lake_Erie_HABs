% this routine combines total CI values for composite images, obtained form
% MERIS data, into one text file

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

datenum_wrapper = @(x)datenum(x,'yyyy-mm-dd');
cc_data=[];
datenums = [];
for dir_ind = 1:length(product_direc)
    
    filename = fullfile(direc,product_direc{dir_ind},'composite_product','cloud_cover.txt');
    fid = fopen(filename,'r');
    data = textscan(fid,'%s%f%s','delimiter','\t','headerlines',1);
    fclose(fid);
    
    cc_data = [cc_data;data{2}];
    
    dates = data{3};
    datenums = [datenums;cellfun(datenum_wrapper,dates)];
end

wfname = 'cloud_cover_combined_MERIS.txt';
filename = fullfile(direc,wfname);
fid = fopen(filename,'w');
fprintf(fid,'%s\t%s\n','Begin_date(yyyy-mm-dd)','cloud_cover(fraction)');
for wind = 1:length(cc_data)
    
    fprintf(fid,'%s\t%d\n',datestr(datenums(wind),'yyyy-mm-dd'),cc_data(wind));
    
end
fclose(fid);


