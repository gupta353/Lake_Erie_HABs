% create a data set for model development

clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/';
agency = 'lake_erie_NOAA';
station = 'erie-cmt_bydate';
met_par = 'daily_wspeed.txt';

% read CI data
fname = 'CI_total_combined_MERIS.txt';
filename = fullfile(direc,'remote_sensing_data',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);

CI_dates = data{1};
CI_vals = data{2};
CI_datenums = cellfun(@(x)datenum(x,'yyyy-mm-dd'),CI_dates);

% read wind-speed data
filename = fullfile(direc,agency,'meteorological_data',station,'daily_wspeed.txt');
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);

ws_dates = data{1};
ws_vals = data{2};
ws_datenums = cellfun(@(x)datenum(x,'dd-mmm-yyyy'),ws_dates);

% read air temperature data
filename = fullfile(direc,agency,'meteorological_data',station,'daily_atemp.txt');
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);

atemp_dates = data{1};
atemp_vals = data{2};
atemp_datenums = cellfun(@(x)datenum(x,'dd-mmm-yyyy'),atemp_dates);

% read total phosphorus data
fname = 'maumee_reconstructed_TP.txt';
filename = fullfile(direc,'HTLP',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f%f%f','delimiter','\t','headerlines',1);
fclose(fid);
TP_dates = data{1};
TP = data{4};
TP_datenums = cellfun(@(x)datenum(x,'mm/dd/yyyy'),TP_dates);

% read total Nitrogen data
fname = 'maumee_reconstructed_TKN_conc.txt';
filename = fullfile(direc,'HTLP',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f%f%f','delimiter','\t','headerlines',1);
fclose(fid);
TKN_dates = data{1};
TKN = data{4};
TKN_datenums = cellfun(@(x)datenum(x,'mm/dd/yyyy'),TP_dates);

% read sreamflow data
fname = 'maumee.txt';
filename = fullfile(direc,'HTLP/streamflow',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);
strm_dates = data{1};
strm_vals = data{2}*0.0283; % conversion to cms
strm_datenums = cellfun(@(x)datenum(x,'mm/dd/yyyy'),strm_dates);

% read secchi depth data
fname = 'secchi_depth_over_pos_CI_total_combined_MERIS.txt';
filename = fullfile(direc,'remote_sensing_data',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);

sd_dates = data{1};
sd_vals = data{2};
sd_datenums = cellfun(@(x)datenum(x,'yyyy-mm-dd'),CI_dates);

%% compute average (or minimum or maximum) of met and TP data over time-period of each CI values (/composite image)
for dind = 1:length(CI_datenums)
    
    ind = find(ws_datenums>=CI_datenums(dind) & ws_datenums<=CI_datenums(dind)+9);
    avg_ws(dind) = min([ws_vals(ind);NaN]);
    ind = find(atemp_datenums>=CI_datenums(dind) & atemp_datenums<=CI_datenums(dind)+9);
    avg_atemp(dind) = max([atemp_vals(ind);NaN]);
    ind = find(TP_datenums>=CI_datenums(dind) & TP_datenums<=CI_datenums(dind)+9);
    avg_TP(dind) = nanmean(TP(ind));
    ind = find(TKN_datenums>=CI_datenums(dind) & TKN_datenums<=CI_datenums(dind)+9);
    avg_TKN(dind) = nanmean(TKN(ind));
    ind = find(strm_datenums>=CI_datenums(dind) & strm_datenums<=CI_datenums(dind)+9);
    avg_strm(dind) = nansum(strm_vals(ind));
    ind = find(sd_datenums>=CI_datenums(dind) & sd_datenums<=CI_datenums(dind)+9);
    avg_sd(dind) = nansum(sd_vals(ind));
    
end
TN_TP_ratio = avg_TKN./avg_TP;

%% compute spring TP and TKN

years = {'2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012'};
years_datenum = cellfun(@(x)datenum(x,'yyyy'),years);

for year_ind = 1:length(years_datenum)-1
    
    
    % total TP load
    begin_date = strcat('01-03-',years{year_ind});
    end_date = strcat('30-06-',years{year_ind});
    begin_datenum = datenum(begin_date,'dd-mm-yyyy');
    end_datenum = datenum(end_date,'dd-mm-yyyy');
    
    ind1 = find(TP_datenums>=begin_datenum & TP_datenums<=end_datenum);
    tot_TP(year_ind) = sum(TP(ind1));
    
    ind1 = find(TKN_datenums>=begin_datenum & TKN_datenums<=end_datenum);
    tot_TKN(year_ind) = sum(TKN(ind1));
end

% add previous 10 years of TP
years = {'2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012'};
years_datenum = cellfun(@(x)datenum(x,'yyyy'),years);

for year_ind = 1:length(years_datenum)-1
    
    
    % total TP load
    begin_year = num2str(str2num(years{year_ind})-10);
    begin_date = strcat('01-01-',begin_year);
    
    end_year = num2str(str2num(years{year_ind})-1);
    end_date = strcat('31-12-',end_year);
    
    begin_datenum = datenum(begin_date,'dd-mm-yyyy');
    end_datenum = datenum(end_date,'dd-mm-yyyy');
    
    ind1 = find(TP_datenums>=begin_datenum & TP_datenums<=end_datenum);
    tot_legacy_TP(year_ind) = sum(TP(ind1));
    
end

%% arrange data in a matrix
data = [CI_datenums,CI_vals,avg_ws',avg_atemp',avg_TP',avg_TKN',TN_TP_ratio',avg_strm',avg_sd'];

%% add previous time-step data
data = sortrows(data);
data = [data,NaN*ones(size(data,1),8)];


for dind = 2:size(data,1)
    
    datenum_tmp = data(dind,1);
    if datenum_tmp - data(dind-1,1)==10
        data(dind,10:end) = data(dind-1,2:9);
    end
    
end

% remove 'seechi depth at current time-step' column
data(:,9) = [];
% add spring TP and TKN
data = [data,zeros(size(data,1),3)];
for dind = 1:size(data,1)
    
    datenum_tmp = data(dind,1);
    year_tmp = datestr(datenum_tmp,'yyyy');
    datenum_tmp = datenum(year_tmp,'yyyy');
    ind = find(years_datenum == datenum_tmp);
    data(dind,end-2) = tot_TP(ind);
    data(dind,end-1) = tot_TKN(ind);
    data(dind,end) = tot_legacy_TP(ind);
    
end


% save data to a textfile
fname = 'model_data.txt';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB','matlab_codes',fname);
fid = fopen(filename,'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','begin_date','CI(t)','min_wind_speed(t)(m/s)','max_air_temperature(t)(\circC)','average_TP(t)(Kg/day)','average_TKN(t)(Kg/day)',...
    'TP_TKN_ratio(t)','average_streamflow(t)(cms)','CI(t-1)','min_wind_speed(t-1)(m/s)','max_air_temperature(t-1)(\circC)','average_TP(t-1)(Kg/day)','average_TKN(t-1)(Kg/day)',...
    'TP_TKN_ratio(t-1)','average_streamflow(t-1)(cms)','secchi_depth(m)','total_spring_TP(Kg)','total_spring_TKN(Kg)','10_year_legacy_TP(Kg)')
for dind = 1:size(data,1)
    
    fprintf(fid,'%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',datestr(data(dind,1),'dd-mmm-yyyy'),data(dind,2:end));

end
fclose(fid);



