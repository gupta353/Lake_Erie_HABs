% this routine combines meteorological data for different years into one
% textfile

clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/lake_erie_NOAA/meteorological_data/erie-cmt_bydate';

fnames = {'wdir_avg_2004.txt','wdir_avg_2005.txt','wdir_avg_2006.txt','wdir_avg_2007.txt','wdir_avg_2008.txt','wdir_avg_2009.txt','wdir_avg_2010.txt',...
    'wdir_avg_2011.txt','wdir_avg_2012.txt','wdir_avg_2013.txt','wdir_avg_2014.txt','wdir_avg_2015.txt'};

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
fname = 'wdir.txt';
filename = fullfile(direc,fname);
fid = fopen(filename,'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n','year','month','day','hour','minute','sec','wind_direction(Deg)');
for wind = 1:size(final_data,1)
    
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%f\n',final_data(wind,:));
    
end
fclose(fid);