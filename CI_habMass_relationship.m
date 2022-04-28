%{
This script explores the relationship between CI and HAB pigments mass per
unit area

Author: Abhinav Gupta (Created: 24 Apr 2022)

%}
clear all
close all
clc

%% read CI data
direc = 'D:/Research/EPA_project/Lake_Erie_HAB/Data/remote_sensing_data/Sentinel';
fname = 'hab_field_sampling_NOAA_locations_CI_value.txt';
filename = fullfile(direc, fname);

fid = fopen(filename, 'r');
data = textscan(fid, '%s%f%f%f%f', 'delimiter', '\t', 'headerlines', 1);
fclose(fid);
dates = data{4};
for i = 1:length(dates)
    datenums_CI(i) = datenum(num2str(dates(i)), 'yyyymmdd');
end
CIs = data{5};
lat_CI = data{2};
lon_CI = data{3};
stn_CI = data{1};

%% read pigment concentration data
direc = 'D:/Research/EPA_project/Lake_Erie_HAB/Data/lake_erie_NOAA/hab_field_sampling_and_nutrient_buoys/1.1/data/0-data';
fname = 'phyco_chl_conc_per_unit_area.txt';
filename = fullfile(direc, fname);
fid = fopen(filename, 'r');
data = textscan(fid, '%s%s%f%f', 'delimiter', '\t', 'headerlines', 1);
fclose(fid);

datenums_hab = cellfun(@(x)datenum(x, 'dd-mmm-yyyy'), data{1});
stn_hab = data{2};
phyco_hab = data{3};
chl_hab = data{4};

%% identify the common dates bewteen two datasets
station_list = {'WE2', 'WE4', 'WE6', 'WE8', 'WE9', 'WE12', 'WE13', 'WE14', 'WE15', 'WE16'};
final_data = [];
for stn_ind = 1:length(station_list)
    stn_tmp = station_list{stn_ind};
    
    ind = find(strcmp(stn_CI ,stn_tmp));
    CIs_tmp = CIs(ind);
    datenums_CI_tmp = datenums_CI(ind);
    
    %
    ind = find(strcmp(stn_hab, stn_tmp));
    datenums_hab_tmp = datenums_hab(ind);
    phyco_hab_tmp = phyco_hab(ind);
    chl_hab_tmp = chl_hab(ind);
    
    for dind = 1:length(datenums_CI_tmp)
        datenum_tmp = datenums_CI_tmp(dind);
        ind = find(datenums_hab_tmp==datenum_tmp);
        if ~isempty(ind)
            final_data = [final_data; [CIs(dind), phyco_hab_tmp(ind), chl_hab_tmp(ind)]];
        end
    end
    scatter(final_data(:,1), final_data(:,3));
end

