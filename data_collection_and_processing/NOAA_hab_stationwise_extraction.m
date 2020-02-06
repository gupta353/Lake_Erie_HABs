% This routine reads the data from
% 'lake_erie_habs_field_sampling_results_2012_2018' and creates a different
% file for each station 

clear all
close all
clc

direc=['D:/Research/EPA_Project/Lake_Erie_HAB/Data/lake_erie_NOAA/'...
    'hab_field_sampling_and_nutrient_buoys/1.1/data/0-data'];

% read data
fname='lake_erie_habs_field_sampling_results_2012_2018.csv';
filename=fullfile(direc,fname);
formatspec=repmat('%s',1,36);
fid=fopen(filename,'r');
data=textscan(fid,formatspec,'delimiter',',');
fclose(fid);

stations=data{2}; stations(1)=[];
station_list=unique(stations);
data=cat(2,data{:});
var_names=data(1,:);
var_names(2)=[];
data(1,:)=[];
data(:,2)=[];


% read data for each station and write into a separate textfile
for station_ind=1:length(station_list)
    
    temp_station=station_list{station_ind};
    ind=find(strcmp(stations,temp_station));
    
    write_fname=strcat('hab_',temp_station,'.txt');
    write_filename=fullfile(direc,write_fname);
    wfid=fopen(write_filename,'wt');
    formatspec=[repmat('%s\t',1,34),'%s\n'];
    fprintf(wfid,formatspec,var_names{1:end});
    
    write_data=data(ind,1:end); write_data=write_data';
    fprintf(wfid,formatspec,write_data{:});
    fclose(wfid);
    
end