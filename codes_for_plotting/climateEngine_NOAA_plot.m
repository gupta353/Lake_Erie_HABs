% this script compares data obtained NOAA-cmt station and climate engine

clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data';
sname = 'atemp_max_NOAA_climate_engine.svg';
axis_text = 'daily max air temperature (\circ C)';
%% read climate engine data
fname = 'ClimateEngine_data_combined.csv';
varName = 'Max_Temperature_degC_';

% read data
filename = fullfile(direc,'climate_engine_data',fname);
data = readtable(filename);
dates = data.DateTime;
varData_climat_engine = data.(varName);
datenums_climat_engine = cellfun(@(x)datenum(x,'mm/dd/yyyy'),dates);

%% read NOAA data
agency = 'lake_erie_NOAA';
station = 'erie-cmt_bydate';
sp_direc = fullfile(direc,agency,'meteorological_data',station);
fname = 'daily_max_atemp.txt';

filename = fullfile(sp_direc,fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);
dates = data{1};
varData_NOAA = data{2};

wrapper = @(x)datenum(x,'dd-mmm-yyyy');
datenums_NOAA = cellfun(wrapper,dates);

%% create a matrix containing NOAA and climate data for each data (replace missing data by NaN)
filled_data = [];
for cind = 1:length(datenums_climat_engine)
    
    datenum_tmp = datenums_climat_engine(cind);
    ind = find(datenums_NOAA==datenum_tmp);
    if ~isempty(ind)
        filled_data(cind,:) = [varData_climat_engine(cind),varData_NOAA(ind)];
    else
        filled_data(cind,:) = [varData_climat_engine(cind),NaN];
    end
        
end

%% plot data
scatter(filled_data(:,1),filled_data(:,2),'filled')
min_val = min(filled_data(:));
max_val = max(filled_data(:));
xlim([min_val max_val]); ylim([min_val max_val]);
hold on;
plot([min_val max_val],[min_val max_val],'linewidth',2,'color','black');
xlabel(['Climate engine ',axis_text],'fontname','arial','fontsize',10);
ylabel(['NOAA cmt ',axis_text],'fontname','arial','fontsize',10);
box('on');
box.linewidth = 2;
set(gca,'fontname','arial','fontsize',10,box);

% save plot
save_filename = fullfile(direc,'climate_engine_data','plots',sname);
fig2svg(save_filename);