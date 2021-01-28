% plot daily scale meteorological data
 clear all
 close all
 clc
 
direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/';
agency = 'lake_erie_NOAA';
station = 'erie-cmt_bydate';

sp_direc = fullfile(direc,agency,'meteorological_data',station);

years = {'2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'};
year_datenum = cellfun(@(x)datenum(x,'yyyy'),years);
%% daily air temperature
%
fname = 'daily_atemp.txt';

filename = fullfile(sp_direc,fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);
dates = data{1};
par = data{2};

wrapper = @(x)datenum(x,'dd-mmm-yyyy');
datenums = cellfun(wrapper,dates);
max_par = max(par);
scatter(datenums,par,5,'filled')
hold on
for y_ind = 1:length(years)
    plot([year_datenum(y_ind) year_datenum(y_ind)],[-30 max_par],'--','linewidth',1,'color','black');
end
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',12,'plotboxaspectratio',[2,1,1],'Xtick',year_datenum,box)
datetick('x','yyyy','keeplimits','keepticks')
xlabel('Year','fontname','arial','fontsize',12)
ylabel('Daily air temperature ( \circ C)','fontname','arial','fontsize',12)
title('Meteorological station: cmt (83.19 W, 41.825 N)','fontname','arial','fontsize',12)

sname = 'daily_air_temperature_cmt.svg';
filename = fullfile(sp_direc,sname);
fig2svg(filename)

clear box
%}
%% wind speed
%{
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
scatter(datenums,par,5,'filled')
hold on
for y_ind = 1:length(years)
    plot([year_datenum(y_ind) year_datenum(y_ind)],[0 max_par],'--','linewidth',2,'color','black');
end
box('on')
box.linewidth=2;
set(gca,'fontname','arial','fontsize',12,'plotboxaspectratio',[2,1,1],'Xtick',year_datenum,box)
datetick('x','yyyy','keeplimits','keepticks')
xlabel('Year','fontname','arial','fontsize',12)
ylabel('Daily wind speed ( m s^{-1})','fontname','arial','fontsize',12)
title('Meteorological station: cmt (83.19 W, 41.825 N)','fontname','arial','fontsize',12)

sname = 'daily_wind-speed_cmt.svg';
filename = fullfile(sp_direc,sname);
fig2svg(filename)

clear box
%}
%% water temperature
%{
fname = 'daily_wtemp.txt';

filename = fullfile(sp_direc,fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);
dates = data{1};
par = data{2};

wrapper = @(x)datenum(x,'dd-mmm-yyyy');
datenums = cellfun(wrapper,dates);
max_par = max(par);
scatter(datenums,par,5,'filled')
hold on
for y_ind = 1:length(years)
    plot([year_datenum(y_ind) year_datenum(y_ind)],[0 max_par],'--','linewidth',2,'color','black');
end
box('on')
box.linewidth = 2;
set(gca,'fontname','arial','fontsize',12,'plotboxaspectratio',[2,1,1],'Xtick',year_datenum,box)
datetick('x','yyyy','keeplimits','keepticks')
xlabel('Year','fontname','arial','fontsize',12)
ylabel('Daily water temperature ( \circ C)','fontname','arial','fontsize',12)
title('Meteorological station: clv (81.69 W, 41.73 N)','fontname','arial','fontsize',12)

sname = 'daily_water_temperature_clv.svg';
filename = fullfile(sp_direc,sname);
fig2svg(filename)

clear box
%}