% This routine plots TKN, TP, and streamflow data at the Maumee station along with TP at LEC stations

clear all
close all
clc

direc_htlp='D:/Research/EPA_Project/Lake_Erie_HAB/Data/HTLP';
save_dir='D:/Research/EPA_Project/Lake_Erie_HAB/Data/plots';

datenum_wrapper=@(x)datenum(x,'mm/dd/yyyy');

% read HTLP data
fname='daily_maumeedata.txt';
filename=fullfile(direc_htlp,fname);
fid=fopen(filename,'r');
formatspec=['%s',repmat('%f',1,6)];
data=textscan(fid,formatspec,'delimiter','\t',...
    'headerlines',1);
fclose(fid);
date_htlp=data{1};
strm_htlp=data{2};
TKN_conc_htlp=data{4};
ind=find(TKN_conc_htlp>0);
date_htlp=date_htlp(ind);
TKN_conc_htlp=TKN_conc_htlp(ind);
strm_htlp=strm_htlp(ind);
TKN_load_htlp=strm_htlp.*TKN_conc_htlp*10^(-3)*24*3600;
datenum_htlp=cellfun(datenum_wrapper,date_htlp);

% plot data
plot(datenum_htlp,strm_htlp);

xlabel('Date','fontname','arial','fontsize',12);
ylabel('streamflow (m^{3} s^{-1})','fontname','arial','fontsize',12);
set(gca,'fontname','arial','fontsize',12,'xlim',[datenum('01-01-1975','dd-mm-yyyy') datenum('31-12-2019','dd-mm-yyyy')],...
    'plotboxaspectratio',[2 1 1])
datetick('x','dd-mmm-yyyy','keeplimits','keepticks');

sname = 'Streamflow_maumee_river.svg';
save_filename=fullfile(save_dir,sname);
fig2svg(save_filename);