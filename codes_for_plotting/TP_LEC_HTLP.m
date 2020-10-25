% This routine plots TP data at the Maumee station along with TP at LEC stations

clear all
close all
clc

direc_lec='D:/Research/EPA_Project/Lake_Erie_HAB/Data/lake_erie_LEC';
direc_htlp='D:/Research/EPA_Project/Lake_Erie_HAB/Data/HTLP';
save_dir='D:/Research/EPA_Project/Lake_Erie_HAB/Data/plots';

LEC_stations={'WQ_4P','WQ_7M','WQ_8M',...
    'WQ_GR1','WQ_MB18','WQ_MB20'};
station_num=6;

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
TP_conc_htlp=data{4};
ind=find(TP_conc_htlp>0);
date_htlp=date_htlp(ind);
TP_conc_htlp=TP_conc_htlp(ind);
strm_htlp=strm_htlp(ind);
TP_load_htlp=strm_htlp.*TP_conc_htlp*10^(-3)*24*3600;
datenum_htlp=cellfun(datenum_wrapper,date_htlp);

% read LEC data
fname_LEC=LEC_stations{station_num};

station_split=strsplit(fname_LEC,'_');
station_name=strcat('LEC,',char(32),station_split{2});

filename=fullfile(direc_lec,fname_LEC);
T=readtable(filename,'delimiter','\t');
date_LEC=T.date;
TP_LEC=T.tp_mg_L_;
datenum_LEC=cellfun(datenum_wrapper,date_LEC);

% plot LEC and HTLP data
%{
[AX,hL,hR]=plotyy(datenum_htlp,TP_load_htlp,datenum_LEC,TP_LEC,@plot,@scatter);
set(hR,'markerfacecolor','r','markeredgecolor','r')
hold on
begin_datenum=datenum('01/01/2000','mm/dd/yyyy'); 
end_datenum=datenum('01/01/2019','mm/dd/yyyy');
set(AX(1),'fontname','arial','fontsize',14,'xlim',[begin_datenum end_datenum]);
set(AX(2),'fontname','arial','fontsize',14,'xAxisLocation','Top',...
    'xlim',[begin_datenum end_datenum],'xTick',[]);
datetick(AX(1),'x','mm/dd/yyyy','keeplimits','keepticks');
datetick(AX(2),'x','mm/dd/yyyy','keeplimits','keepticks');
% datetick(AX(2),'x','mm/dd/yyyy','keeplimits','keepticks');
xlabel(AX(1),'Date','fontname','arial','fontsize',14);
ylabel(AX(1),'TP, average daily load at Mammee inlet (Kg day^{-1})','fontname','arial','fontsize',14);
ylabel(AX(2),'TP (mg L^{-1})','fontname','arial','fontsize',14);
title(station_name,'fontname','arial','fontsize',14);
xticklabel_rotate([],15);

save_filename=fullfile(save_dir,strcat(LEC_stations{station_num},'.svg'));
fig2svg(save_filename);
%}

% plot LEC data
scatter(datenum_LEC,TP_LEC,'filled')
begin_datenum=datenum('01/01/2000','mm/dd/yyyy'); 
end_datenum=datenum('01/01/2019','mm/dd/yyyy');
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',12,'xlim',[begin_datenum end_datenum],box);
xlabel('Date','fontname','arial','fontsize',12);
datetick('x','dd-mmm-yyyy','keeplimits','keepticks');
ylabel('Total phopshorus (mg L^{-1})','fontname','arial','fontsize',12);
title(['Station: ',station_name],'fontname','arial','fontsize',12);

sname = ['TP_',station_name,'.svg'];
filename = fullfile(save_dir,sname);
fig2svg(filename);