% plots of total CI against meteorological data

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

% read total phosphorus data
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
strm_vals = data{2};
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
    
    %         ind = find(ws_datenums>=CI_datenums(dind) & ws_datenums<=CI_datenums(dind)+9);
    %         avg_ws(dind) = min([ws_vals(ind);NaN]);
    %         ind = find(atemp_datenums>=CI_datenums(dind) & atemp_datenums<=CI_datenums(dind)+9);
    %         avg_atemp(dind) = max([atemp_vals(ind);NaN]);
    ind = find(TP_datenums>=CI_datenums(dind) & TP_datenums<=CI_datenums(dind)+9);
    avg_TP(dind) = nanmean(TP(ind));
    ind = find(TKN_datenums>=CI_datenums(dind) & TKN_datenums<=CI_datenums(dind)+9);
    avg_TKN(dind) = nanmean(TKN(ind));
    %     ind = find(strm_datenums>=CI_datenums(dind) & strm_datenums<=CI_datenums(dind)+9);
    %     avg_strm(dind) = nansum(strm_vals(ind));
    
end
TN_TP = avg_TKN./avg_TP;
% scatter plot
%{
scatter(avg_met,CI_vals,5,'filled')
xlabel('10-day average air temperature (\circC)','fontname','arial','fontsize',12)
ylabel('10-day composite CI total','fontname','arial','fontsize',12)
box('on')
box.linewidth=2;
set(gca,'fontname','arial','fontsize',12,box);
title('Meteorological station: cmt (81.69 W, 41.83 N)','fontname','arial','fontsize',12)
break
sname = 'CI_vs_atemp_cmt.svg';
filename = fullfile(direc,'remote_sensing_data',sname);
fig2svg(filename)
%}

%%
%{
years = {'2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012'};
years_datenum = cellfun(@(x)datenum(x,'yyyy'),years);
months = {'05','06','07','08','09','10'};


for year_ind = 1:length(years_datenum)-1
    
    ind1 = find(CI_datenums>=years_datenum(year_ind)& CI_datenums<years_datenum(year_ind+1));
    for dind = 1:length(months)
        datenums(dind) = datenum(str2double(years(year_ind)),str2double(months(dind)),01);
    end
    
    [Ax,L1,L2]=plotyy(CI_datenums(ind1),avg_strm(ind1),CI_datenums(ind1),log(CI_vals(ind1)));
    set(L1,'linestyle','--','marker','o','markerfacecolor','b','markersize',5,'linewidth',2)
    set(L2,'marker','o','markerfacecolor','g','markersize',5,'linewidth',2)
    box('on')
    box.linewidth=2;
    set(Ax(2),'xAxislocation','top','Xtick',[],'fontname','arial','fontsize',12)
    set(Ax(1),'Xtick',datenums,'fontname','arial','fontsize',12,box)
    ylabel(Ax(1),'Minimum wind speed (m s^{-1})','fontname','arial','fontsize',12)
    ylabel(Ax(2),'log(CI)','fontname','arial','fontsize',12)
    xlabel(Ax(1),'Date','fontname','arial','fontsize',12)
    title(['year = ',years{year_ind}],'fontname','arial','fontsize',12)
    
    datetick(Ax(1),'x','dd-mmm','keeplimits','keepticks')
    datetick(Ax(2),'x','dd-mmm','keeplimits','keepticks')
    clear box
    
    %
    %         sname = ['min_wind_',years{year_ind},'.svg'];
    %         filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB','matlab_codes/plots',sname);
    %     %     print(filename,'-r300','-djpeg');
    %         fig2svg(filename);
    close all
end
%}
%% 3D plots (wind speed, temperature, and CI)
%{
plot_data = [CI_vals,avg_atemp',avg_ws'];
col = repmat('b',size(plot_data,1),1);
bubbleplot(plot_data(:,2),plot_data(:,3),plot_data(:,1),plot_data(:,1),col,[],'MarkerSizeLimits', [5 100]);
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',12,box)
grid on
xlabel('Air temperature (\circC)','fontname','arial','fontsize',12);
ylabel('Wind speed (m s^{-1})','fontname','arial','fontsize',12);
zlabel('CI','fontname','arial','fontsize',12)
%}

%% compute time to peak with 1st May of the year as reference and compute total TP load in Lake from Jan to April
%
years = {'2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012'};
years_datenum = cellfun(@(x)datenum(x,'yyyy'),years);
months = {'05','06','07','08','09','10'};

for year_ind = 1:length(years_datenum)-1
    
    % time to peak
    ind1 = find(CI_datenums>=years_datenum(year_ind)& CI_datenums<years_datenum(year_ind+1));
    comp_datenums = (CI_datenums(ind1(1)):10:CI_datenums(ind1(end)))';
    CI_tmp = NaN*ones(length(comp_datenums),1);
    
    for ii = 1:length(ind1)
        ind = ind1(ii);
        CI_tmp(comp_datenums==CI_datenums(ind)) = CI_vals(ind);
    end
    
    begin_date = strcat('01-05-',years{year_ind});
    begin_datenum = datenum(begin_date,'dd-mm-yyyy');
    mean_CI(year_ind) = nanmean(CI_tmp);
    max_CI(year_ind) = max(CI_tmp);
    max_CI_str{year_ind} = num2str(round(max_CI(year_ind)));
    tpeak(year_ind) = comp_datenums(CI_tmp==max_CI(year_ind));
    tpeak(year_ind) = tpeak(year_ind) - begin_datenum;
    int_CI(year_ind) = CI_tmp(1);       % initial CI
    int_max_CI(year_ind) = max(CI_tmp(1:3));       % initial max CI
    int_CI_str{year_ind} = num2str(round(CI_tmp(1)));       % initial CI string
    
 %   total TP load
    begin_date = strcat('01-01-',years{year_ind});
    end_date = strcat('30-06-',years{year_ind});
    begin_datenum = datenum(begin_date,'dd-mm-yyyy');
    end_datenum = datenum(end_date,'dd-mm-yyyy');
    
    ind1 = find(TP_datenums>=begin_datenum & TP_datenums<=end_datenum);
    tot_TP(year_ind) = sum(TP(ind1));
    tot_TP_str{year_ind} = num2str(round(tot_TP(year_ind)));
    
    ind1 = find(TKN_datenums>=begin_datenum & TKN_datenums<=end_datenum);
    tot_TKN(year_ind) = sum(TKN(ind1));
end

% plot of time to annnual-peak-of-CI vs. met factors
bubbleplot(tot_TP,max_CI);

xlabel('Total phosphorus from Jan to Jun','fontname','arial','fontsize',12);
ylabel('Maximum annual CI','fontname','arial','fontsize',12);
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',12,box);

sname = 'max_annual_CI_vs_TP_Jan_Jun.svg';
filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB','matlab_codes/plots',sname);
fig2svg(filename);
%}

%% plot of CI vs. secchi depth
%{
years = {'2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012'};
years_datenum = cellfun(@(x)datenum(x,'yyyy'),years);
months = {'05','06','07','08','09','10'};
xtick_labels = cellfun(@(x)datenum(x,'mm'),months);

for year_ind = 1:length(years_datenum)-1
    
    ind1 = find(CI_datenums>=years_datenum(year_ind)& CI_datenums<years_datenum(year_ind+1));
    ind2 = find(sd_datenums>=years_datenum(year_ind)& sd_datenums<years_datenum(year_ind+1));
    
    % time-series plot
%     [Ax,L1,L2]=plotyy(sd_datenums(ind2),sd_vals(ind2),CI_datenums(ind1),CI_vals(ind1));
%     set(L1,'linestyle','--','marker','o','markerfacecolor','b','markersize',5,'linewidth',2)
%     set(L2,'marker','o','markerfacecolor','g','markersize',5,'linewidth',2)
%     box('on')
%     box.linewidth=2;
%     set(Ax(2),'xAxislocation','top','Xtick',[],'fontname','arial','fontsize',12)
%     set(Ax(1),'Xticklabel',xtick_labels,'fontname','arial','fontsize',12,box)
%     ylabel(Ax(2),'CI','fontname','arial','fontsize',12)
%     ylabel(Ax(1),'Average secchi depth (m)','fontname','arial','fontsize',12)
%     xlabel(Ax(1),'Date','fontname','arial','fontsize',12)
%     title(['year = ',years{year_ind}],'fontname','arial','fontsize',12)
%     
%     datetick(Ax(1),'x','dd-mmm','keeplimits','keepticks')
%     datetick(Ax(2),'x','dd-mmm','keeplimits','keepticks')
%     clear box
%     pause;
    %
    
    % scatter plots
%         [int_datenums] = intersect(CI_datenums(ind1),sd_datenums(ind2));
%         plot_data = [];
%         for dind = 1:length(int_datenums)
%             iCI = find(CI_datenums==int_datenums(dind));
%             isd = find(sd_datenums==int_datenums(dind));
%             plot_data = [plot_data;[sd_vals(isd),CI_vals(iCI)]];
%         end
%         scatter(plot_data(1:end-2,1),plot_data(3:end,2))
%         pause;
%         clear plot_data
    
    %
%     sname = ['sd__over_pos_CI_',years{year_ind},'.svg'];
%     filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB','matlab_codes/plots',sname);
%     print(filename,'-r300','-djpeg');
%     fig2svg(filename);
%     close all
    
    
end
%}