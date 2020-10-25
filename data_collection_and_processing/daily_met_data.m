% this routine computes daily wind speed, daily average water and air
% temperatures

clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/';
agency = 'lake_erie_NOAA';
station = 'erie-cmt_bydate';

sp_direc = fullfile(direc,agency,'meteorological_data',station);

fname = 'wspped_avg.txt';
wfname = 'daily_wspped_avg.txt';
header='wind_spped(m/s)';

filename = fullfile(sp_direc,fname);
fid = fopen(filename,'r');
data = textscan(fid,'%f%f%f%f%f%f%f','delimiter','\t','headerlines',1);
fclose(fid);
data = cat(2,data{:});

% remove observations of 1929
% ind = find(data(:,1)==1929);
% data(ind,:) = [];

% compute datenums
datenums = datenum(data(:,1),data(:,2),data(:,3));
uni_datenums = unique(datenums,'stable');
uni_datenums(isnan(uni_datenums))=[];

% remove all the uni_datenums that do not start with 7
datenum_tmp = floor(uni_datenums/10^5);
uni_datenums(datenum_tmp~=7)=[];

% remove all the data after labeled to be collected 2020
ind = find(uni_datenums>datenum('2020-01-01','yyyy-mm-dd'));
uni_datenums(ind)=[];

for date_ind = 1:length(uni_datenums)
    
    ind = find(datenums==uni_datenums(date_ind));
    tmp_data = data(ind,7);
    tmp_data(isnan(tmp_data)) = [];
    avg_par(date_ind) = mean(tmp_data);

end

% remove any negative values from the averaged data
ind = find(avg_par<0);
avg_par(ind) = [];
uni_datenums(ind) = [];

% write data to a text file
wfilename = fullfile(sp_direc,wfname);
wfid = fopen(wfilename,'w');
fprintf(wfid,'%s\t%s\n','Date(yyyy-mm-dd)',header);

for wind = 1:length(avg_par)
    
    fprintf(wfid,'%s\t%f\n',datestr(uni_datenums(wind)),avg_par(wind));
    
end
fclose(wfid);
