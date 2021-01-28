% plot daily scale meteorological data year wise
clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/';
agency = 'lake_erie_NOAA';
station = 'erie-cmt_bydate';

sp_direc = fullfile(direc,agency,'meteorological_data',station);

years = {'2002','2003','2004','2005','2006','2007','2008','2009','2010','2011'};
year_datenum = cellfun(@(x)datenum(x,'yyyy'),years);

fname = 'daily_wspeed.txt';

filename = fullfile(sp_direc,fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);
dates = data{1};
par = data{2};

wrapper = @(x)datenum(x,'dd-mmm-yyyy');
datenums = cellfun(wrapper,dates);
max_par = max(par);

count=0;
for year_ind = 1:length(years)
    begin_datenum = datenum(['05/01/',years{year_ind}],'mm/dd/yyyy');
    end_datenum = datenum(['05/30/',years{year_ind}],'mm/dd/yyyy');
    ind = find(datenums>=begin_datenum & datenums<=end_datenum);
    scatter(datenums(ind)-begin_datenum+1,par(ind),50,'filled');
    hold on
    
    if ~isempty(ind)
        count = count+1;
        lgd_labels{count} = years{year_ind};
    end
end
legend(lgd_labels,'fontname','arial','fontsize',12,'Location','Northwest','Orientation','horizontal');
legend('boxoff');
box('on')
box.linewidth = 2;
xlabel('Day in the month of May ','fontname','arial','fontsize',12)
ylabel('Daily wind speed (m s^{-1})','fontname','arial','fontsize',12)
title('Meteorological station: cmt (813.193 W, 41.825 N)','fontname','arial','fontsize',12)

sname = 'daily_wspeed_May_cmt.svg';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB','matlab_codes/plots',sname);
fig2svg(filename);