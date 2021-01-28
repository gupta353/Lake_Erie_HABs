% this routine combines meteorological data for different years into one
% textfile

clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/lake_erie_NOAA/meteorological_data/erie-cmt_bydate/fine_timescale_data';

fnames = {'atemp_2004.txt','atemp_2005.txt','atemp_2006.txt','atemp_2007.txt','atemp_2008.txt','atemp_2009.txt','atemp_2010.txt',...
    'atemp_2011.txt','atemp_2012.txt','atemp_2013.txt','atemp_2014.txt','atemp_2015.txt','atemp_2016.txt',...
    'atemp_2017.txt','atemp_2018.txt','atemp_2019.txt','atemp_2020.txt'};

final_data=[];
for fname_ind = 1:length(fnames)
    
    fname = fnames{fname_ind};
    filename = fullfile(direc,fname);
    fid = fopen(filename,'r');
    data = textscan(fid,'%f%f%f%f%f%f%f','delimiter','\t','headerlines',1);
    fclose(fid);
    data = cat(2,data{:});
    final_data = [final_data;data];
end

% write data to a textfile
fname = 'atemp.txt';
filename = fullfile(direc,fname);
fid = fopen(filename,'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n','year','month','day','hour','minute','sec','atemp(DegC)');
for wind = 1:size(final_data,1)
    
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%f\n',final_data(wind,:));
    
end
fclose(fid);