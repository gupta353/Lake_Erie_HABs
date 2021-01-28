% this script processes data doanloaded from NOAA GLERL, for a particular
% year
% Note: This script is written to process data downlaoded in the sensor format

clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/lake_erie_NOAA/meteorological_data/erie-cmt_bydate';
year  = '2016';
subdir = 'metcsi';

% meteorological variable to read and its column number in the raw files
varName = 'wdir_avg';
col = 13; % metadata file says that the column number is actually 'col-1'; however, an additional comma is introduced after the date while downloading the data

% list all the files in the relevant folder
f_direc = fullfile(direc,year,subdir);
list = dir(f_direc);

write_data = [];
for list_ind = 1:size(list,1)
    
    % read the files in for loop
    fname = list(list_ind).name;
    
    if ~strcmp(fname,{'.','..'})
        filename = fullfile(f_direc,fname);
        fid = fopen(filename,'r');
        formatspec = repmat('%s',1,24);
        data = textscan(fid,formatspec,'headerlines',1,'delimiter',',');
        fclose(fid);
        
        % read dates corresponding to each timestamp
        dates = data{1};
        wrapper_dates = @(x)strsplit(x,{'"','-'});
        dates = cellfun(wrapper_dates,dates,'UniformOutput',false);
        dates = cat(1,dates{:});
        dates(:,1) = [];
        dates = cellfun(@str2num,dates);
        
        % read time (HH:MM:SS) corresponding to each timestamp
        time_hhmmss = data{2};
        wrapper_time = @(x)strsplit(x,{'"',':'});
        time_hhmmss = cellfun(wrapper_time,time_hhmmss,'UniformOutput',false);
        time_hhmmss = cat(1,time_hhmmss{:});
        time_hhmmss(:,1) = [];
        time_hhmmss(:,4) = [];
        time_hhmmss = cellfun(@str2num,time_hhmmss);
        
        % read variable data
        varData = data{col};
        varData = cellfun(@str2num,varData);
        
        write_data = [write_data;[dates,time_hhmmss,varData]];
    end
end

% write data to a textfile
fname = [varName,'_',year,'.txt'];
filename = fullfile(direc,'fine_timescale_data',fname);
fid = fopen(filename,'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n','year','month','day','hour','minute','sec','wspeed_avg(m/s)');
for wind = 1:size(write_data,1)
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%d\n',num2str(write_data(wind,1)),num2str(write_data(wind,2)),num2str(write_data(wind,3)),num2str(write_data(wind,4)),num2str(write_data(wind,5)),num2str(write_data(wind,6)),write_data(wind,7));
end
fclose(fid);
