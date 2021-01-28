% this script processes data doanloaded from NOAA GLERL, for a particular
% year
% Note: This script is written to process data downlaoded in the csv format

clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/lake_erie_NOAA/meteorological_data';
year = '2020';
varName = 'wdir_rmy';
wfname = 'wdir_avg';
numHeaderLines = 29;

% read data
fname = [varName,'_',year,'.csv'];
filename = fullfile(direc,'erie-cmt_csv',year,fname);
fid = fopen(filename,'r');
data=textscan(fid,'%s%f%f','headerlines',numHeaderLines,'delimiter',',');
fclose(fid);

dates = data{1};
varData = data{2};

wrapper = @(x)strsplit(x,{'/',char(32),':'});
dates = cellfun(wrapper,dates,'UniformOutput',false);
dates = cat(1,dates{:});

% write data to a textfile
fname = [wfname,'_',year,'.txt'];
filename = fullfile(direc,'erie-cmt_bydate/fine_timescale_data',fname);
fid = fopen(filename,'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n','year','month','day','hour','minute','sec','wind_direction(Deg)');
for wind = 1:size(dates,1)
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%f\n',dates{wind,1},dates{wind,2},dates{wind,3},dates{wind,4},dates{wind,5},'0',varData(wind));
end
fclose(fid);