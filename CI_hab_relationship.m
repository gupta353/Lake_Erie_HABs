%{
This script explores the relationship between CI and HAB pigment concentrations

Author: Abhinav Gupta (Created: 22 Apr 2022)

%}
clear all
close all
clc

% read CI data
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

%% read microcystin concentration data
direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/lake_erie_NOAA/hab_field_sampling_and_nutrient_buoys/1.1/data/0-data';
fname = 'lake_erie_habs_field_sampling_results_2012_2018.txt';
filename = fullfile(direc, fname);
fid = fopen(filename, 'r');
formatspec = repmat('%s',1,36);
data = textscan(fid, formatspec, 'delimiter', {'','\t'}, 'headerlines', 1);
fclose(fid);
dates = data{1};
datenums_hab = cellfun(@(x)datenum(x, 'mm/dd/yyyy'), dates);
stn_hab = data{2};
type = data{5};
phyco_str = data{23};
chl_str = data{24};
depth_str = data{4};

% convert strings to float
for ind = 1:length(phyco_str)
    
    if strcmp(phyco_str{ind}, '')
        phyco(ind) = NaN;
    elseif strcmp(phyco_str{ind}(1), '<')
       phyco(ind) = NaN;
    else
        phyco(ind) = str2num(phyco_str{ind});
    end
    
    if strcmp(chl_str{ind}, '')
        chl(ind) = NaN;
    elseif strcmp(chl_str{ind}(1), '<')
       chl(ind) = NaN;
    else
        chl(ind) = str2num(chl_str{ind});
    end
    
    if strcmp(depth_str{ind}, '')
        depth(ind) = NaN;
    elseif strcmp(depth_str{ind}(1), '<')
       depth(ind) = NaN;
    else
        depth(ind) = str2num(depth_str{ind});
    end
    
end

% kepp only surface samples
ind = find(strcmp(type, 'Surface'));
datenums_hab = datenums_hab(ind);
stn_hab = stn_hab(ind);
phyco = phyco(ind);
chl = chl(ind);
depth = depth(ind);
%% extract data for each station in a for loop from CI and hab files
station_list = {'WE2', 'WE4', 'WE6', 'WE8', 'WE9', 'WE12', 'WE13', 'WE14', 'WE15', 'WE16'};

plot_data = [];
for sind = 1:length(station_list)
    
    stn = station_list{sind};
    
    % extract CI data
    ind = find(strcmp(stn_CI, stn));
    CIs_tmp = CIs(ind);
    datenums_tmp_CI = datenums_CI(ind);
    
    % extract hab data
    ind = find(strcmp(stn_hab, stn));
    datenums_tmp_hab = datenums_hab(ind);
    phyco_tmp = phyco(ind);
    chl_tmp = chl(ind);
    depth_tmp = depth(ind);
    
    for dind = 1:length(datenums_tmp_CI)
        ind = find(datenums_tmp_hab == datenums_tmp_CI(dind));
        if ~isempty(ind)
            plot_data = [plot_data; [CIs_tmp(dind), phyco_tmp(ind), chl_tmp(ind)]];
        end
    end
    
end

col_ind = 3;
scatter(plot_data(:,1), plot_data(:,col_ind), 'filled')
set(gca, 'xscale', 'linear')
hold on

% fit linear model
ind = find(~isnan(plot_data(:,1)));
mdl = fitlm(plot_data(ind,1), plot_data(ind,col_ind));
slope = mdl.Coefficients.Estimate(2);
inter = mdl.Coefficients.Estimate(1);
plot(plot_data(ind,1), inter + slope*plot_data(ind,1), 'black');
xlabel('CI');
ylabel(['Phycocyanin concentration (', char(181), 'g L^{-1})'])
box('on')
box.linewidth = 1;
set(gca, 'fontname', 'Arial', 'fontsize', 12, 'GridLineStyle', '--', box);
grid on

% save plot
direc_save = 'D:/Research/EPA_project/Lake_Erie_HAB/matlab_codes/plots_08_28_2021/CI_hab_relationship';
sname = 'CI_chl_surface.svg';
filename = fullfile(direc_save, sname);
saveas(gcf, filename);

% save model
sname = 'CI_chl_surface.mat';
filename = fullfile(direc_save, sname);
save(filename, 'mdl')