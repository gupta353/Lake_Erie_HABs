% create a data set for model development

clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/';
agency = 'lake_erie_NOAA';
station = 'erie-cmt_bydate';
met_par = 'daily_wspeed.txt';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read CI data
fname = 'CI_total_MERIS_SEN.txt';
filename = fullfile(direc,'remote_sensing_data',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);

CI_dates = data{1};
CI_vals = data{2};
CI_datenums = cellfun(@(x)datenum(x,'yyyy-mm-dd'),CI_dates);

% remove teh data corresponding to year 2020
ind = find(CI_datenums>=datenum('01/01/2020','mm/dd/yyyy'));
CI_datenums(ind) = [];
CI_vals(ind) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read wind-speed data
filename = fullfile(direc,agency,'meteorological_data',station,'daily_wspeed_reconstructed.txt');
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);

ws_dates = data{1};
ws_vals = data{2};
ws_datenums = cellfun(@(x)datenum(x,'mm/dd/yyyy'),ws_dates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read mean air temperature data
filename = fullfile(direc,agency,'meteorological_data',station,'daily_mean_atemp_reconstructed.txt');
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);

atemp_dates = data{1};
mean_atemp_vals = data{2};
mean_atemp_datenums = cellfun(@(x)datenum(x,'mm/dd/yyyy'),atemp_dates);

% read maximum air temperature data
filename = fullfile(direc,agency,'meteorological_data',station,'daily_max_atemp_reconstructed.txt');
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);

atemp_dates = data{1};
max_atemp_vals = data{2};
max_atemp_datenums = cellfun(@(x)datenum(x,'mm/dd/yyyy'),atemp_dates);

% read minimum air temperature data
filename = fullfile(direc,agency,'meteorological_data',station,'daily_min_atemp_reconstructed.txt');
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);

atemp_dates = data{1};
min_atemp_vals = data{2};
min_atemp_datenums = cellfun(@(x)datenum(x,'mm/dd/yyyy'),atemp_dates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read total phosphorus data from Maumee river
fname = 'maumee_reconstructed_TP_conc.txt';
filename = fullfile(direc,'HTLP_updated',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f%f%f','delimiter','\t','headerlines',1);
fclose(fid);
TP_dates = data{1};
TP_maumee = data{4};
TP_datenums_maumee = cellfun(@(x)datenum(x,'mm/dd/yyyy'),TP_dates);

% read total phosphorus data from Raisin river
fname = 'raisin_reconstructed_TP_conc.txt';
filename = fullfile(direc,'HTLP_updated',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f%f%f','delimiter','\t','headerlines',1);
fclose(fid);
TP_dates = data{1};
TP_raisin = data{4};
TP_datenums_raisin = cellfun(@(x)datenum(x,'mm/dd/yyyy'),TP_dates);

% read total phosphorus data from sandusky river
fname = 'sandusky_reconstructed_TP_conc.txt';
filename = fullfile(direc,'HTLP_updated',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f%f%f','delimiter','\t','headerlines',1);
fclose(fid);
TP_dates = data{1};
TP_sandusky = data{4};
TP_datenums_sandusky = cellfun(@(x)datenum(x,'mm/dd/yyyy'),TP_dates);

% read total phosphorus data from cuyahoga river
fname = 'cuyahoga_reconstructed_TP_conc.txt';
filename = fullfile(direc,'HTLP_updated',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f%f%f','delimiter','\t','headerlines',1);
fclose(fid);
TP_dates = data{1};
TP_cuyahoga = data{4};
TP_datenums_cuyahoga = cellfun(@(x)datenum(x,'mm/dd/yyyy'),TP_dates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read total Nitrogen data from maumee river
fname = 'maumee_reconstructed_TKN_conc.txt';
filename = fullfile(direc,'HTLP_updated',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f%f%f','delimiter','\t','headerlines',1);
fclose(fid);
TKN_dates = data{1};
TKN_maumee = data{4};
TKN_datenums_maumee = cellfun(@(x)datenum(x,'mm/dd/yyyy'),TP_dates);

% read total Nitrogen data from raisin river
fname = 'raisin_reconstructed_TKN_conc.txt';
filename = fullfile(direc,'HTLP_updated',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f%f%f','delimiter','\t','headerlines',1);
fclose(fid);
TKN_dates = data{1};
TKN_raisin = data{4};
TKN_datenums_raisin = cellfun(@(x)datenum(x,'mm/dd/yyyy'),TP_dates);

% read total Nitrogen data from sandusky river
fname = 'sandusky_reconstructed_TKN_conc.txt';
filename = fullfile(direc,'HTLP_updated',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f%f%f','delimiter','\t','headerlines',1);
fclose(fid);
TKN_dates = data{1};
TKN_sandusky = data{4};
TKN_datenums_sandusky = cellfun(@(x)datenum(x,'mm/dd/yyyy'),TP_dates);

% read total Nitrogen data from cuyahoga river
fname = 'cuyahoga_reconstructed_TKN_conc.txt';
filename = fullfile(direc,'HTLP_updated',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f%f%f','delimiter','\t','headerlines',1);
fclose(fid);
TKN_dates = data{1};
TKN_cuyahoga = data{4};
TKN_datenums_cuyahoga = cellfun(@(x)datenum(x,'mm/dd/yyyy'),TP_dates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read sreamflow data from maumee river
fname = 'maumee.txt';
filename = fullfile(direc,'HTLP_updated/streamflow',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);
strm_dates = data{1};
strm_vals_maumee = data{2}*0.0283; % conversion to cms
strm_datenums_maumee = cellfun(@(x)datenum(x,'mm/dd/yyyy'),strm_dates);

% read sreamflow data from raisin river
fname = 'raisin.txt';
filename = fullfile(direc,'HTLP_updated/streamflow',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);
strm_dates = data{1};
strm_vals_raisin = data{2}*0.0283; % conversion to cms
strm_datenums_raisin = cellfun(@(x)datenum(x,'mm/dd/yyyy'),strm_dates);

% read sreamflow data from sandusky river
fname = 'sandusky.txt';
filename = fullfile(direc,'HTLP_updated/streamflow',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);
strm_dates = data{1};
strm_vals_sandusky = data{2}*0.0283; % conversion to cms
strm_datenums_sandusky = cellfun(@(x)datenum(x,'mm/dd/yyyy'),strm_dates);

% read sreamflow data from cuyahoga river
fname = 'cuyahoga.txt';
filename = fullfile(direc,'HTLP_updated/streamflow',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);
strm_dates = data{1};
strm_vals_cuyahoga = data{2}*0.0283; % conversion to cms
strm_datenums_cuyahoga = cellfun(@(x)datenum(x,'mm/dd/yyyy'),strm_dates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read water level data
fname = 'water_level_mean.txt';
filename = fullfile(direc,'water_level',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);

WL_dates = data{1};
WL_vals = data{2};
WL_datenums = cellfun(@(x)datenum(x,'mm/dd/yyyy'),WL_dates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read solar radiation data
fname = 'solar_radiation.csv';
filename = fullfile(direc,'climate_engine_data',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter',',','headerlines',1);
fclose(fid);

SR_dates = data{1};
SR_vals = data{2};
SR_datenums = cellfun(@(x)datenum(x,'mm/dd/yyyy'),SR_dates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read secchi depth data
fname = 'secchi_depth_over_pos_CI_total_MERIS_SEN.txt';
filename = fullfile(direc,'remote_sensing_data',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);

sd_dates = data{1};
sd_vals = data{2};
sd_datenums = cellfun(@(x)datenum(x,'yyyy-mm-dd'),CI_dates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute average (or minimum or maximum) of met and TP data over time-period of each CI values (/composite image)
for dind = 1:length(CI_datenums)
    
    % wind speed
    ind = find(ws_datenums>=CI_datenums(dind) & ws_datenums<=CI_datenums(dind)+9);
    min_ws(dind) = min([ws_vals(ind);NaN]);
    max_ws(dind) = max([ws_vals(ind);NaN]);
    avg_ws(dind) = nanmean([ws_vals(ind);NaN]);
    
    % max air temperature
    ind = find(max_atemp_datenums>=CI_datenums(dind) & max_atemp_datenums<=CI_datenums(dind)+9);
    max_atemp(dind) = max([max_atemp_vals(ind);NaN]);
    
    % min air temperature
    ind = find(min_atemp_datenums>=CI_datenums(dind) & min_atemp_datenums<=CI_datenums(dind)+9);
    min_atemp(dind) = min([min_atemp_vals(ind);NaN]);
    
    % mean air temperature
    ind = find(mean_atemp_datenums>=CI_datenums(dind) & mean_atemp_datenums<=CI_datenums(dind)+9);
    avg_atemp(dind) = nanmean([mean_atemp_vals(ind);NaN]);
    
    % solar radiation
    ind = find(SR_datenums>=CI_datenums(dind) & SR_datenums<=CI_datenums(dind)+9);
    avg_SR(dind) = nanmean([SR_vals(ind);NaN]);
    
    % water level
    ind = find(WL_datenums>=CI_datenums(dind) & WL_datenums<=CI_datenums(dind)+9);
    avg_WL(dind) = nanmean([WL_vals(ind);NaN]);
    %% average TP within 10 days
    % from maumee river
    ind = find(TP_datenums_maumee>=CI_datenums(dind) & TP_datenums_maumee<=CI_datenums(dind)+9);
    avg_TP_maumee(dind) = nanmean(TP_maumee(ind));
    
    % from raisin river
    ind = find(TP_datenums_raisin>=CI_datenums(dind) & TP_datenums_raisin<=CI_datenums(dind)+9);
    avg_TP_raisin(dind) = nanmean(TP_raisin(ind));
    
    % from sandusky river
    ind = find(TP_datenums_sandusky>=CI_datenums(dind) & TP_datenums_sandusky<=CI_datenums(dind)+9);
    avg_TP_sandusky(dind) = nanmean(TP_sandusky(ind));
    
    % from cuyahoga river
    ind = find(TP_datenums_cuyahoga>=CI_datenums(dind) & TP_datenums_cuyahoga<=CI_datenums(dind)+9);
    avg_TP_cuyahoga(dind) = nanmean(TP_cuyahoga(ind));
    
    %% average TKN within 10 days
    % from maumee river
    ind = find(TKN_datenums_maumee>=CI_datenums(dind) & TKN_datenums_maumee<=CI_datenums(dind)+9);
    avg_TKN_maumee(dind) = nanmean(TKN_maumee(ind));
    
    % from raisin river
    ind = find(TKN_datenums_raisin>=CI_datenums(dind) & TKN_datenums_raisin<=CI_datenums(dind)+9);
    avg_TKN_raisin(dind) = nanmean(TKN_raisin(ind));
    
    % from sandusky river
    ind = find(TKN_datenums_sandusky>=CI_datenums(dind) & TKN_datenums_sandusky<=CI_datenums(dind)+9);
    avg_TKN_sandusky(dind) = nanmean(TKN_sandusky(ind));
    
    % from cuyahoga river
    ind = find(TKN_datenums_cuyahoga>=CI_datenums(dind) & TKN_datenums_cuyahoga<=CI_datenums(dind)+9);
    avg_TKN_cuyahoga(dind) = nanmean(TKN_cuyahoga(ind));
    
    %% total streamflow within 10 days
    % from maumee river
    ind = find(strm_datenums_maumee>=CI_datenums(dind) & strm_datenums_maumee<=CI_datenums(dind)+9);
    avg_strm_maumee(dind) = nansum(strm_vals_maumee(ind));
    
    % from raisin river
    ind = find(strm_datenums_raisin>=CI_datenums(dind) & strm_datenums_raisin<=CI_datenums(dind)+9);
    avg_strm_raisin(dind) = nansum(strm_vals_raisin(ind));
    
    % from sandusky river
    ind = find(strm_datenums_sandusky>=CI_datenums(dind) & strm_datenums_sandusky<=CI_datenums(dind)+9);
    avg_strm_sandusky(dind) = nansum(strm_vals_sandusky(ind));
    
    % from cuyahoga river
    ind = find(strm_datenums_cuyahoga>=CI_datenums(dind) & strm_datenums_cuyahoga<=CI_datenums(dind)+9);
    avg_strm_cuyahoga(dind) = nansum(strm_vals_cuyahoga(ind));
    
    %% average Secchi depth 10 days
    ind = find(sd_datenums>=CI_datenums(dind) & sd_datenums<=CI_datenums(dind)+9);
    avg_sd(dind) = nanmean(sd_vals(ind));   
    
end
TN_TP_ratio_maumee = avg_TKN_maumee./avg_TP_maumee;
TN_TP_ratio_raisin = avg_TKN_raisin./avg_TP_raisin;
TN_TP_ratio_sandusky = avg_TKN_sandusky./avg_TP_sandusky;
TN_TP_ratio_cuyahoga = avg_TKN_sandusky./avg_TP_cuyahoga;

%% compute spring TP and TKN
%
years = {'2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2016','2017','2018','2019'};
years_datenum = cellfun(@(x)datenum(x,'yyyy'),years);

for year_ind = 1:length(years_datenum)
    
    
    % total TP load
    begin_date = strcat('01-03-',years{year_ind});
    end_date = strcat('30-06-',years{year_ind});
    begin_datenum = datenum(begin_date,'dd-mm-yyyy');
    end_datenum = datenum(end_date,'dd-mm-yyyy');
    
    % from maumee river
    ind1 = find(TP_datenums_maumee>=begin_datenum & TP_datenums_maumee<=end_datenum);
    tot_TP_maumee(year_ind) = sum(TP_maumee(ind1));
    ind1 = find(TKN_datenums_maumee>=begin_datenum & TKN_datenums_maumee<=end_datenum);
    tot_TKN_maumee(year_ind) = sum(TKN_maumee(ind1));
    
    % from raisin river
    ind1 = find(TP_datenums_raisin>=begin_datenum & TP_datenums_raisin<=end_datenum);
    tot_TP_raisin(year_ind) = sum(TP_raisin(ind1));
    ind1 = find(TKN_datenums_raisin>=begin_datenum & TKN_datenums_raisin<=end_datenum);
    tot_TKN_raisin(year_ind) = sum(TKN_raisin(ind1));
    
    % from sandusky river
    ind1 = find(TP_datenums_sandusky>=begin_datenum & TP_datenums_sandusky<=end_datenum);
    tot_TP_sandusky(year_ind) = sum(TP_sandusky(ind1));
    ind1 = find(TKN_datenums_sandusky>=begin_datenum & TKN_datenums_sandusky<=end_datenum);
    tot_TKN_sandusky(year_ind) = sum(TKN_sandusky(ind1));
    
    % from cuyahoga river
    ind1 = find(TP_datenums_cuyahoga>=begin_datenum & TP_datenums_cuyahoga<=end_datenum);
    tot_TP_cuyahoga(year_ind) = sum(TP_cuyahoga(ind1));
    ind1 = find(TKN_datenums_cuyahoga>=begin_datenum & TKN_datenums_cuyahoga<=end_datenum);
    tot_TKN_cuyahoga(year_ind) = sum(TKN_cuyahoga(ind1));
    
end
%}

%% compute TP and TKN within last 120 days
%
for dind = 1:length(CI_datenums)
    
    % total TP load
    begin_datenum = CI_datenums(dind)-120;
    end_datenum = CI_datenums(dind)-1;
    
    % from maumee river
    ind1 = find(TP_datenums_maumee>=begin_datenum & TP_datenums_maumee<=end_datenum);
    tot_TP_120_maumee(dind) = sum(TP_maumee(ind1));
    ind1 = find(TKN_datenums_maumee>=begin_datenum & TKN_datenums_maumee<=end_datenum);
    tot_TKN_120_maumee(dind) = sum(TKN_maumee(ind1));
    
    % from raisin river
    ind1 = find(TP_datenums_raisin>=begin_datenum & TP_datenums_raisin<=end_datenum);
    tot_TP_120_raisin(dind) = sum(TP_raisin(ind1));
    ind1 = find(TKN_datenums_raisin>=begin_datenum & TKN_datenums_raisin<=end_datenum);
    tot_TKN_120_raisin(dind) = sum(TKN_raisin(ind1));
    
    % from sandusky river
    ind1 = find(TP_datenums_sandusky>=begin_datenum & TP_datenums_sandusky<=end_datenum);
    tot_TP_120_sandusky(dind) = sum(TP_sandusky(ind1));
    ind1 = find(TKN_datenums_sandusky>=begin_datenum & TKN_datenums_sandusky<=end_datenum);
    tot_TKN_120_sandusky(dind) = sum(TKN_sandusky(ind1));
    
    % from cuyahoga river
    ind1 = find(TP_datenums_cuyahoga>=begin_datenum & TP_datenums_cuyahoga<=end_datenum);
    tot_TP_120_cuyahoga(dind) = sum(TP_cuyahoga(ind1));
    ind1 = find(TKN_datenums_cuyahoga>=begin_datenum & TKN_datenums_cuyahoga<=end_datenum);
    tot_TKN_120_cuyahoga(dind) = sum(TKN_cuyahoga(ind1));
    
end

%% add previous 10 years of TP
years = {'2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2016','2017','2018','2019'};
years_datenum = cellfun(@(x)datenum(x,'yyyy'),years);

for year_ind = 1:length(years_datenum)
    
    
    % total TP load
    begin_year = num2str(str2num(years{year_ind})-10);
    begin_date = strcat('01-01-',begin_year);
    
    end_year = num2str(str2num(years{year_ind})-1);
    end_date = strcat('31-12-',end_year);
    
    begin_datenum = datenum(begin_date,'dd-mm-yyyy');
    end_datenum = datenum(end_date,'dd-mm-yyyy');
    
    % from maumee river
    ind1 = find(TP_datenums_maumee>=begin_datenum & TP_datenums_maumee<=end_datenum);
    tot_legacy_TP_maumee(year_ind) = sum(TP_maumee(ind1));
    
    % from raisin river
    ind1 = find(TP_datenums_raisin>=begin_datenum & TP_datenums_raisin<=end_datenum);
    tot_legacy_TP_raisin(year_ind) = sum(TP_raisin(ind1));
    
    % from sandusky river
    ind1 = find(TP_datenums_sandusky>=begin_datenum & TP_datenums_sandusky<=end_datenum);
    tot_legacy_TP_sandusky(year_ind) = sum(TP_sandusky(ind1));
    
    % from cuyahoga river
    ind1 = find(TP_datenums_cuyahoga>=begin_datenum & TP_datenums_cuyahoga<=end_datenum);
    tot_legacy_TP_cuyahoga(year_ind) = sum(TP_cuyahoga(ind1));
    
end

%% compute TP, TKN and streamflow integrated in time within last 20,30,40,50,and 60 days
%{
for dind = 1:length(CI_datenums)
    
    %% average TP within last 20,30, 40,50,and 60 days
    % from maumee river
    ind = find(TP_datenums_maumee<CI_datenums(dind) & TP_datenums_maumee>=CI_datenums(dind)-30);
    avg_TP_maumee_30(dind) = nanmean(TP_maumee(ind));
    
    ind = find(TP_datenums_maumee<CI_datenums(dind) & TP_datenums_maumee>=CI_datenums(dind)-40);
    avg_TP_maumee_40(dind) = nanmean(TP_maumee(ind));
    
    ind = find(TP_datenums_maumee<CI_datenums(dind) & TP_datenums_maumee>=CI_datenums(dind)-50);
    avg_TP_maumee_50(dind) = nanmean(TP_maumee(ind));
    
    ind = find(TP_datenums_maumee<CI_datenums(dind) & TP_datenums_maumee>=CI_datenums(dind)-60);
    avg_TP_maumee_60(dind) = nanmean(TP_maumee(ind));
    
    % from raisin river
    ind = find(TP_datenums_raisin<CI_datenums(dind) & TP_datenums_raisin>=CI_datenums(dind)-30);
    avg_TP_raisin_30(dind) = nanmean(TP_raisin(ind));
    
    ind = find(TP_datenums_raisin<CI_datenums(dind) & TP_datenums_raisin>=CI_datenums(dind)-40);
    avg_TP_raisin_40(dind) = nanmean(TP_raisin(ind));
    
    ind = find(TP_datenums_raisin<CI_datenums(dind) & TP_datenums_raisin>=CI_datenums(dind)-50);
    avg_TP_raisin_50(dind) = nanmean(TP_raisin(ind));
    
    ind = find(TP_datenums_raisin<CI_datenums(dind) & TP_datenums_raisin>=CI_datenums(dind)-60);
    avg_TP_raisin_60(dind) = nanmean(TP_raisin(ind));
    
    % from sandusky river
    ind = find(TP_datenums_sandusky<CI_datenums(dind) & TP_datenums_sandusky>=CI_datenums(dind)-30);
    avg_TP_sandusky_30(dind) = nanmean(TP_sandusky(ind));
    
    ind = find(TP_datenums_sandusky<CI_datenums(dind) & TP_datenums_sandusky>=CI_datenums(dind)-40);
    avg_TP_sandusky_40(dind) = nanmean(TP_sandusky(ind));
    
    ind = find(TP_datenums_sandusky<CI_datenums(dind) & TP_datenums_sandusky>=CI_datenums(dind)-50);
    avg_TP_sandusky_50(dind) = nanmean(TP_sandusky(ind));
    
    ind = find(TP_datenums_sandusky<CI_datenums(dind) & TP_datenums_sandusky>=CI_datenums(dind)-60);
    avg_TP_sandusky_60(dind) = nanmean(TP_sandusky(ind));
    
    % from cuyahoga river
    ind = find(TP_datenums_cuyahoga<CI_datenums(dind) & TP_datenums_cuyahoga>=CI_datenums(dind)-30);
    avg_TP_cuyahoga_30(dind) = nanmean(TP_cuyahoga(ind));
    
    ind = find(TP_datenums_cuyahoga<CI_datenums(dind) & TP_datenums_cuyahoga>=CI_datenums(dind)-40);
    avg_TP_cuyahoga_40(dind) = nanmean(TP_cuyahoga(ind));
    
    ind = find(TP_datenums_cuyahoga<CI_datenums(dind) & TP_datenums_cuyahoga>=CI_datenums(dind)-50);
    avg_TP_cuyahoga_50(dind) = nanmean(TP_cuyahoga(ind));
    
    ind = find(TP_datenums_cuyahoga<CI_datenums(dind) & TP_datenums_cuyahoga>=CI_datenums(dind)-60);
    avg_TP_cuyahoga_60(dind) = nanmean(TP_cuyahoga(ind));
    
     %% average TKN within last 20,30, 40,50,and 60 days
    % from maumee river
    ind = find(TKN_datenums_maumee<CI_datenums(dind) & TKN_datenums_maumee>=CI_datenums(dind)-30);
    avg_TKN_maumee_30(dind) = nanmean(TKN_maumee(ind));
    
    ind = find(TKN_datenums_maumee<CI_datenums(dind) & TKN_datenums_maumee>=CI_datenums(dind)-40);
    avg_TKN_maumee_40(dind) = nanmean(TKN_maumee(ind));
    
    ind = find(TKN_datenums_maumee<CI_datenums(dind) & TKN_datenums_maumee>=CI_datenums(dind)-50);
    avg_TKN_maumee_50(dind) = nanmean(TKN_maumee(ind));
    
    ind = find(TKN_datenums_maumee<CI_datenums(dind) & TKN_datenums_maumee>=CI_datenums(dind)-60);
    avg_TKN_maumee_60(dind) = nanmean(TKN_maumee(ind));
    
    % from raisin river
    ind = find(TKN_datenums_raisin<CI_datenums(dind) & TKN_datenums_raisin>=CI_datenums(dind)-30);
    avg_TKN_raisin_30(dind) = nanmean(TKN_raisin(ind));
    
    ind = find(TKN_datenums_raisin<CI_datenums(dind) & TKN_datenums_raisin>=CI_datenums(dind)-40);
    avg_TKN_raisin_40(dind) = nanmean(TKN_raisin(ind));
    
    ind = find(TKN_datenums_raisin<CI_datenums(dind) & TKN_datenums_raisin>=CI_datenums(dind)-50);
    avg_TKN_raisin_50(dind) = nanmean(TKN_raisin(ind));
    
    ind = find(TKN_datenums_raisin<CI_datenums(dind) & TKN_datenums_raisin>=CI_datenums(dind)-60);
    avg_TKN_raisin_60(dind) = nanmean(TKN_raisin(ind));
    
    % from sandusky river
    ind = find(TKN_datenums_sandusky<CI_datenums(dind) & TKN_datenums_sandusky>=CI_datenums(dind)-30);
    avg_TKN_sandusky_30(dind) = nanmean(TKN_sandusky(ind));
    
    ind = find(TKN_datenums_sandusky<CI_datenums(dind) & TKN_datenums_sandusky>=CI_datenums(dind)-40);
    avg_TKN_sandusky_40(dind) = nanmean(TKN_sandusky(ind));
    
    ind = find(TKN_datenums_sandusky<CI_datenums(dind) & TKN_datenums_sandusky>=CI_datenums(dind)-50);
    avg_TKN_sandusky_50(dind) = nanmean(TKN_sandusky(ind));
    
    ind = find(TKN_datenums_sandusky<CI_datenums(dind) & TKN_datenums_sandusky>=CI_datenums(dind)-60);
    avg_TKN_sandusky_60(dind) = nanmean(TKN_sandusky(ind));
    
    % from cuyahoga river
    ind = find(TKN_datenums_cuyahoga<CI_datenums(dind) & TKN_datenums_cuyahoga>=CI_datenums(dind)-30);
    avg_TKN_cuyahoga_30(dind) = nanmean(TKN_cuyahoga(ind));
    
    ind = find(TKN_datenums_cuyahoga<CI_datenums(dind) & TKN_datenums_cuyahoga>=CI_datenums(dind)-40);
    avg_TKN_cuyahoga_40(dind) = nanmean(TKN_cuyahoga(ind));
    
    ind = find(TKN_datenums_cuyahoga<CI_datenums(dind) & TKN_datenums_cuyahoga>=CI_datenums(dind)-50);
    avg_TKN_cuyahoga_50(dind) = nanmean(TKN_cuyahoga(ind));
    
    ind = find(TKN_datenums_cuyahoga<CI_datenums(dind) & TKN_datenums_cuyahoga>=CI_datenums(dind)-60);
    avg_TKN_cuyahoga_60(dind) = nanmean(TKN_cuyahoga(ind));
    
    %% average streamflows within last 20,30, 40,50,and 60 days
    % from maumee river
    ind = find(strm_datenums_maumee<CI_datenums(dind) & strm_datenums_maumee>=CI_datenums(dind)-30);
    avg_strm_maumee_30(dind) = nanmean(strm_vals_maumee(ind));
    
    ind = find(strm_datenums_maumee<CI_datenums(dind) & strm_datenums_maumee>=CI_datenums(dind)-40);
    avg_strm_maumee_40(dind) = nanmean(strm_vals_maumee(ind));
    
    ind = find(strm_datenums_maumee<CI_datenums(dind) & strm_datenums_maumee>=CI_datenums(dind)-50);
    avg_strm_maumee_50(dind) = nanmean(strm_vals_maumee(ind));
    
    ind = find(strm_datenums_maumee<CI_datenums(dind) & strm_datenums_maumee>=CI_datenums(dind)-60);
    avg_strm_maumee_60(dind) = nanmean(strm_vals_maumee(ind));
    
    % from raisin river
    ind = find(strm_datenums_raisin<CI_datenums(dind) & strm_datenums_raisin>=CI_datenums(dind)-30);
    avg_strm_raisin_30(dind) = nanmean(strm_vals_raisin(ind));
    
    ind = find(strm_datenums_raisin<CI_datenums(dind) & strm_datenums_raisin>=CI_datenums(dind)-40);
    avg_strm_raisin_40(dind) = nanmean(strm_vals_raisin(ind));
    
    ind = find(strm_datenums_raisin<CI_datenums(dind) & strm_datenums_raisin>=CI_datenums(dind)-50);
    avg_strm_raisin_50(dind) = nanmean(strm_vals_raisin(ind));
    
    ind = find(strm_datenums_raisin<CI_datenums(dind) & strm_datenums_raisin>=CI_datenums(dind)-60);
    avg_strm_raisin_60(dind) = nanmean(strm_vals_raisin(ind));
    
    % from sandusky river
    ind = find(strm_datenums_sandusky<CI_datenums(dind) & strm_datenums_sandusky>=CI_datenums(dind)-30);
    avg_strm_sandusky_30(dind) = nanmean(strm_vals_sandusky(ind));
    
    ind = find(strm_datenums_sandusky<CI_datenums(dind) & strm_datenums_sandusky>=CI_datenums(dind)-40);
    avg_strm_sandusky_40(dind) = nanmean(strm_vals_sandusky(ind));
    
    ind = find(strm_datenums_sandusky<CI_datenums(dind) & strm_datenums_sandusky>=CI_datenums(dind)-50);
    avg_strm_sandusky_50(dind) = nanmean(strm_vals_sandusky(ind));
    
    ind = find(strm_datenums_sandusky<CI_datenums(dind) & strm_datenums_sandusky>=CI_datenums(dind)-60);
    avg_strm_sandusky_60(dind) = nanmean(strm_vals_sandusky(ind));
    
    % from cuyahoga river
    ind = find(strm_datenums_cuyahoga<CI_datenums(dind) & strm_datenums_cuyahoga>=CI_datenums(dind)-30);
    avg_strm_cuyahoga_30(dind) = nanmean(strm_vals_cuyahoga(ind));
    
    ind = find(strm_datenums_cuyahoga<CI_datenums(dind) & strm_datenums_cuyahoga>=CI_datenums(dind)-40);
    avg_strm_cuyahoga_40(dind) = nanmean(strm_vals_cuyahoga(ind));
    
    ind = find(strm_datenums_cuyahoga<CI_datenums(dind) & strm_datenums_cuyahoga>=CI_datenums(dind)-50);
    avg_strm_cuyahoga_50(dind) = nanmean(strm_vals_cuyahoga(ind));
    
    ind = find(strm_datenums_cuyahoga<CI_datenums(dind) & strm_datenums_cuyahoga>=CI_datenums(dind)-60);
    avg_strm_cuyahoga_60(dind) = nanmean(strm_vals_cuyahoga(ind));
end
%}
%% arrange data in a matrix
%
data = [CI_datenums,CI_vals,min_ws',avg_ws',max_ws',min_atemp',avg_atemp',max_atemp',...
    avg_TP_maumee',avg_TP_raisin',avg_TP_sandusky',avg_TP_cuyahoga',...
    avg_TKN_maumee',avg_TKN_raisin',avg_TKN_sandusky',avg_TKN_cuyahoga',...
    avg_strm_maumee',avg_strm_raisin',avg_strm_sandusky',avg_strm_cuyahoga',...
    TN_TP_ratio_maumee',TN_TP_ratio_raisin',TN_TP_ratio_sandusky',TN_TP_ratio_cuyahoga',avg_SR',avg_WL',avg_sd'];

%% add previous time-step data
data = sortrows(data);                      % arrange the data in increasing order of dates
data = [data,NaN*ones(size(data,1),26)];     

for dind = 2:size(data,1)
    
    datenum_tmp = data(dind,1);
    if datenum_tmp - data(dind-1,1)==10
        data(dind,28:end) = data(dind-1,2:27);
    end
    
end

%% remove 'seechi depth at current time-step' column
data(:,27) = [];

%% add spring TP and TKN
data = [data,zeros(size(data,1),12)];
for dind = 1:size(data,1)
    
    datenum_tmp = data(dind,1);
    year_tmp = datestr(datenum_tmp,'yyyy');
    datenum_tmp = datenum(year_tmp,'yyyy');
    ind = find(years_datenum == datenum_tmp);

    data(dind,end-11:end) = [tot_TP_maumee(ind),tot_TP_raisin(ind),tot_TP_sandusky(ind),tot_TP_cuyahoga(ind),...
        tot_TKN_maumee(ind),tot_TKN_raisin(ind),tot_TKN_sandusky(ind),tot_TKN_cuyahoga(ind),...
        tot_legacy_TP_maumee(ind),tot_legacy_TP_raisin(ind),tot_legacy_TP_sandusky(ind),tot_legacy_TP_cuyahoga(ind)];
    
end

%% compute average (or minimum or maximum) of met and TP data which have maximum correlation with CI values
for dind = 1:length(CI_datenums)
       
    % max air temperature
    ind = find(max_atemp_datenums>=CI_datenums(dind)-30 & max_atemp_datenums<=CI_datenums(dind)-21);
    max_atemp_corr(dind) = max([max_atemp_vals(ind);NaN]);
    
    % min air temperature
    ind = find(min_atemp_datenums>=CI_datenums(dind)-30 & min_atemp_datenums<=CI_datenums(dind)-21);
    min_atemp_corr(dind) = min([min_atemp_vals(ind);NaN]);
    
    % mean air temperature
    ind = find(mean_atemp_datenums>=CI_datenums(dind)-30 & mean_atemp_datenums<=CI_datenums(dind)-21);
    avg_atemp_corr(dind) = nanmean([mean_atemp_vals(ind);NaN]);
    
    % solar radiation
    ind = find(SR_datenums>=CI_datenums(dind)-60 & SR_datenums<=CI_datenums(dind)-30);
    avg_SR_corr(dind) = nanmean([SR_vals(ind);NaN]);
    
    % water level
    ind = find(WL_datenums>=CI_datenums(dind)-60 & WL_datenums<=CI_datenums(dind)-50);
    avg_WL_corr(dind) = nanmean([WL_vals(ind);NaN]);
    
    % Average TP from maumee river 
    ind = find(TP_datenums_maumee>=CI_datenums(dind)-30 & TP_datenums_maumee<=CI_datenums(dind)-1);
    avg_TP_maumee_30(dind) = nanmean(TP_maumee(ind));
    
    % Average TP from raisin river
    ind = find(TP_datenums_raisin>=CI_datenums(dind)-30 & TP_datenums_raisin<=CI_datenums(dind)-1);
    avg_TP_raisin_30(dind) = nanmean(TP_raisin(ind));
    
    % Average TKN from maumee river
    ind = find(TKN_datenums_maumee>=CI_datenums(dind)-30 & TKN_datenums_maumee<=CI_datenums(dind)-1);
    avg_TKN_maumee_30(dind) = nanmean(TKN_maumee(ind));
    
    % Average TKN from raisin river
    ind = find(TKN_datenums_raisin>=CI_datenums(dind)-30 & TKN_datenums_raisin<=CI_datenums(dind)-1);
    avg_TKN_raisin_30(dind) = nanmean(TKN_raisin(ind));
    
end
data = [data,max_atemp_corr',min_atemp_corr',avg_atemp_corr',avg_SR_corr',avg_WL_corr',avg_TP_maumee_30',avg_TP_raisin_30',avg_TKN_maumee_30',avg_TKN_raisin_30'];

%% Add two time-step lag CI values to dataset
data = [data,NaN*ones(size(data,1),1)];     

for dind = 3:size(data,1)
    
    datenum_tmp = data(dind,1);
    if datenum_tmp - data(dind-2,1)==20
        data(dind,74) = data(dind-2,2);
    end
    
end

%% add previous 120 days total TP and TKN to dataset
data = [data,tot_TP_120_maumee',tot_TP_120_raisin',tot_TP_120_sandusky',tot_TP_120_cuyahoga',tot_TKN_120_maumee',tot_TKN_120_raisin',tot_TKN_120_sandusky',tot_TKN_120_cuyahoga'];

%% Corresponding time-step of the year (1 may to 10 may is 1st time-step, 11th may to 20th may is second time-step etc)
for dind = 1:length(CI_datenums)
    year = datestr(CI_datenums(dind),'yyyy');
    first_datenum = datenum(['01-May-',year],'dd-mmm-yyyy');
    time_step(dind) = (CI_datenums(dind)-first_datenum)/10+1;
end
data = [data,time_step'];
%% save data to a textfile

fname = 'model_data_10_corr_CI_lag_included_and_TP_TKN_previous_120_and_time_step_included.txt';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB','matlab_codes',fname);
fid = fopen(filename,'w');
fprintf(fid,[repmat('%s\t',1,82),'%s\n'],'begin_date','CI(t)','min_wind_speed(t)(m/s)','avg_wind_speed(t)(m/s)','max_wind_speed(t)(m/s)',...
    'min_air_temperature(t)(\circC)','avg_air_temperature(t)(\circC)','max_air_temperature(t)(\circC)',...
    'avg_TP_maumee(t)(Kg/day)','avg_TP_raisin(t)(Kg/day)','avg_TP_sandusky(t)(Kg/day)','avg_TP_cuyahoga(t)(Kg/day)',...
    'avg_TKN_maumee(t)(Kg/day)','avg_TKN_raisin(t)(Kg/day)','avg_TKN_sandusky(t)(Kg/day)','avg_TKN_cuyahoga(t)(Kg/day)',...
    'avg_streamflow_maumee(t)(cms)','avg_streamflow_raisin(t)(cms)','avg_streamflow_sandusky(t)(cms)','avg_streamflow_cuyahoga(t)(cms)',...
    'TP_TKN_ratio_maumee(t)','TP_TKN_ratio_raisin(t)','TP_TKN_ratio_sandusky(t)','TP_TKN_ratio_cuyahoga(t)',...
    'avg_solar_radiation(t)(W/m2)','avg_water_level(t)(m)',...
    'CI(t-1)','min_wind_speed(t-1)(m/s)','avg_wind_speed(t-1)(m/s)','max_wind_speed(t-1)(m/s)',...
    'min_air_temperature(t-1)(\circC)','avg_air_temperature(t-1)(\circC)','max_air_temperature(t-1)(\circC)',...
    'avg_TP_maumee(t-1)(Kg/day)','avg_TP_raisin(t-1)(Kg/day)','avg_TP_sandusky(t-1)(Kg/day)','avg_TP_cuyahoga(t-1)(Kg/day)',...
    'avg_TKN_maumee(t-1)(Kg/day)','avg_TKN_raisin(t-1)(Kg/day)','avg_TKN_sandusky(t-1)(Kg/day)','avg_TKN_cuyahoga(t-1)(Kg/day)',...
    'average_streamflow_maumee(t-1)(cms)','average_streamflow_raisin(t-1)(cms)','average_streamflow_sandusky(t-1)(cms)','average_streamflow_cuyahoga(t-1)(cms)',...
    'TP_TKN_ratio_maumee(t-1)','TP_TKN_ratio_raisin(t-1)','TP_TKN_ratio_sandusky(t-1)','TP_TKN_ratio_cuyahoga(t-1)',...
    'avg_solar_radiation(t-1)(W/m2)','avg_water_level(t-1)(m)','secchi_depth(t-1)(m)',...
    'spring_TP_maumee(Kg)','spring_TP_raisin(Kg)','spring_TP_sandusky(Kg)','spring_TP_cuyahoga(Kg)',...
    'spring_TKN_maumee(Kg)','spring_TKN_raisin(Kg)','spring_TKN_sandusky(Kg)','spring_TKN_cuyahoga(Kg)',...
    'Legacy_TP_maumee(Kg)','Legacy_TP_raisin(Kg)','Legacy_TP_sandusky(Kg)','Legacy_TP_cuyahoga(Kg)','max_atemp_lag_3(\circC)','min_atemp_lag_3(\circC)',...
    'avg_atemp_lag_3(\circC)','solar_radiation_30_lag_3(W/m^2)','avg_water_level_lag_6(m)','avg_TP_maumee_30(Kg/day)','avg_TP_raisin_30(Kg/day)',...
    'avg_TKN_maumee_30(Kg/day)','avg_TKN_raisin_30(Kg/day)','CI(t-2)','tot_TP_120_maumee(Kg)','tot_TP_120_raisin(Kg)','tot_TP_120_sandusky(Kg)','tot_TP_120_cuyahoga(Kg)',...
    'tot_TKN_120_maumee(Kg)','tot_TKN_120_raisin(Kg)','tot_TKN_120_sandusky(Kg)','tot_TKN_120_cuyahoga(Kg)','time_step_of_the_year');

for dind = 1:size(data,1)
    
    fprintf(fid,['%s\t',repmat('%f\t',1,81),'%f\n'],datestr(data(dind,1),'dd-mmm-yyyy'),data(dind,2:end));
    
end
fclose(fid);
%}

%% arrange time-integrated data into a matrix and write to a textfile
%{
data = [CI_vals,avg_TP_maumee_30',avg_TP_maumee_40',avg_TP_maumee_50',avg_TP_maumee_60',avg_TP_raisin_30',avg_TP_raisin_40',avg_TP_raisin_50',...
    avg_TP_raisin_60',avg_TP_sandusky_30',avg_TP_sandusky_40',avg_TP_sandusky_50',avg_TP_sandusky_60',avg_TP_cuyahoga_30',avg_TP_cuyahoga_40',...
    avg_TP_cuyahoga_50',avg_TP_cuyahoga_60',avg_TKN_maumee_30',avg_TKN_maumee_40',avg_TKN_maumee_50',avg_TKN_maumee_60',avg_TKN_raisin_30',...
    avg_TKN_raisin_40',avg_TKN_raisin_50',avg_TKN_raisin_60',avg_TKN_sandusky_30',avg_TKN_sandusky_40',avg_TKN_sandusky_50',avg_TKN_sandusky_60',...
    avg_TKN_cuyahoga_30',avg_TKN_cuyahoga_40',avg_TKN_cuyahoga_50',avg_TKN_cuyahoga_60',avg_strm_maumee_30',avg_strm_maumee_40',avg_strm_maumee_50',...
    avg_strm_maumee_60',avg_strm_raisin_30',avg_strm_raisin_40',avg_strm_raisin_50',avg_strm_raisin_60',avg_strm_sandusky_30',avg_strm_sandusky_40',...
    avg_strm_sandusky_50',avg_strm_sandusky_60',avg_strm_cuyahoga_30',avg_strm_cuyahoga_40',avg_strm_cuyahoga_50',avg_strm_cuyahoga_60'];

% write data to a textfile
fname = 'CI_time_integrated_data.txt';
filename = fullfile(direc,fname);
fid = fopen(filename,'w');
formatspec = [repmat('%s\t',1,48),'%s\n'];
fprintf(fid,formatspec,'CI_vals','avg_TP_maumee_30','avg_TP_maumee_40','avg_TP_maumee_50','avg_TP_maumee_60','avg_TP_raisin_30','avg_TP_raisin_40','avg_TP_raisin_50',...
    'avg_TP_raisin_60','avg_TP_sandusky_30','avg_TP_sandusky_40','avg_TP_sandusky_50','avg_TP_sandusky_60','avg_TP_cuyahoga_30','avg_TP_cuyahoga_40',...
    'avg_TP_cuyahoga_50','avg_TP_cuyahoga_60','avg_TKN_maumee_30','avg_TKN_maumee_40','avg_TKN_maumee_50','avg_TKN_maumee_60','avg_TKN_raisin_30',...
    'avg_TKN_raisin_40','avg_TKN_raisin_50','avg_TKN_raisin_60','avg_TKN_sandusky_30','avg_TKN_sandusky_40','avg_TKN_sandusky_50','avg_TKN_sandusky_60',...
    'avg_TKN_cuyahoga_30','avg_TKN_cuyahoga_40','avg_TKN_cuyahoga_50','avg_TKN_cuyahoga_60','avg_strm_maumee_30','avg_strm_maumee_40','avg_strm_maumee_50',...
    'avg_strm_maumee_60','avg_strm_raisin_30','avg_strm_raisin_40','avg_strm_raisin_50','avg_strm_raisin_60','avg_strm_sandusky_30','avg_strm_sandusky_40',...
    'avg_strm_sandusky_50','avg_strm_sandusky_60','avg_strm_cuyahoga_30','avg_strm_cuyahoga_40','avg_strm_cuyahoga_50','avg_strm_cuyahoga_60');

formatspec = [repmat('%f\t',1,48),'%f\n'];
for data_ind = 1:size(data,1)
    
    fprintf(fid,formatspec,data(data_ind,:));
    
end
fclose(fid);
%}