% this script plots climate engine data

clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/climate_engine_data';
fname = 'ClimateEngine_data_combined.csv';
varName = 'Max_Temperature_degC_';
ylabel_text = 'Maximum daily temperature (\circ C)';
sname = 'atemp_max.svg';

% read data
filename = fullfile(direc,fname);
data = readtable(filename);
dates = data.DateTime;
varData = data.(varName);

% convert date to datenum
datenums = cellfun(@(x)datenum(x,'mm/dd/yyyy'),dates);
yearnums = datenum(num2str([2002:2020]'),'yyyy');

% plot data
plot(datenums,varData);
hold on
for yind = 1:length(yearnums)
    plot([yearnums(yind) yearnums(yind)],[0 max(varData)],'--','color','black')
end

set(gca,'xlim',[min(datenums) max(datenums)],'ylim',[0 max(varData)],'fontname','arial','fontsize',10,'plotboxaspectratio',[2 1 1])
datetick('x','yyyy','keepticks')
xlabel('Day','fontname','arial','fontsize',10);
ylabel(ylabel_text,'fontname','arial','fontsize',10);

% save plot
save_filename = fullfile(direc,'plots',sname);
fig2svg(save_filename);