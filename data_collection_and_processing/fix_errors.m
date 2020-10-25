% this routine some of the error in the files 'atemp_[year]','wdir_[year]',
% and 'wspped_avg_[year]',
% 

clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/lake_erie_NOAA/meteorological_data/erie-cmt_bydate';
fname = 'wspped_avg_2010.txt';

%% remove NaNs from meteorological data
%{
filename = fullfile(direc,fname);
fid = fopen(filename,'r');
header = fgetl(fid);
data = textscan(fid,'%f%f%f%f%f%f%f','delimiter','\t');
fclose(fid);
data = cat(2,data{:});

ind = find(isnan(data(:,3)));

for find = 1:length(ind)
    
    tmp_ind = ind(find);
    data(tmp_ind,1:3) = data(tmp_ind-1,1:3);
end

% write data to a textfile
filename = fullfile(direc,fname);       % file with the same name as the original file
wfid = fopen(filename,'w');
fprintf(wfid,'%s\n',header);
for find = 1:size(data,1)
    fprintf(wfid,'%d\t%d\t%d\t%d\t%d\t%d\t%d\n',data(find,:));
end
fclose(wfid);
%}

%% remove outliers from windspeed data