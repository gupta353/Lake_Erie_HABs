% Reconstruction of maumee river nutient load data using RVM

clear all
close all
clc

direc_htlp='D:/Research/EPA_Project/Lake_Erie_HAB/Data/HTLP';
datenum_wrapper=@(x)datenum(x,'mm/dd/yyyy');
strsplit_wrapper=@(x)strsplit(x,'/');

% read wq concentration data
fname_htlp='daily_maumeedata.txt';
filename=fullfile(direc_htlp,fname_htlp);
data=readtable(filename,'delimiter','\t');
date_wq=data.(1);
datenum_wq=cellfun(datenum_wrapper,date_wq);
wq_conc=data.TP_Mg_LAsP;
wq_conc(wq_conc<0)=NaN;

% read streamflow data
fname='maumee.txt';
filename=fullfile(direc_htlp,'streamflow',fname);
data=readtable(filename,'delimiter','\t');
date_strm=data.(1);
datenum_strm=cellfun(datenum_wrapper,date_strm);
strm=data.(2);              % streamflow in cfs

% find common datenums in wq and strm data
[~,i_wq,i_strm]=intersect(datenum_wq,datenum_strm);
wq_conc=wq_conc(i_wq);
datenum_wq=datenum_wq(i_wq);
date_wq=date_wq(i_wq);

% fill the missing dates with NaNs in wq data and pick up streamflow data
% corresponding towq data
datenum_range=datenum_wq(1):datenum_wq(end);
[~,iwq,~]=intersect(datenum_range,datenum_wq);
wq_conc_tmp=nan(length(datenum_range),1);
wq_conc_tmp(iwq)=wq_conc;
[~,~,istrm]=intersect(datenum_range,datenum_strm);
strm=strm(istrm);

% split the date into yyyy, mm, and dd
for date_ind=1:length(datenum_range)
    
    date{date_ind}=datestr(datenum_range(date_ind),'mm/dd/yyyy');
    date_wq_split_tmp=strsplit(date{date_ind},'/');
    date_wq_split(date_ind,:)=[str2double(date_wq_split_tmp(3)),...
        str2double(date_wq_split_tmp(1)),str2double(date_wq_split_tmp(2))];
end

% reconstruct wq data
modelNo=3;  % which LAODEST model to use
unitno='1'; % 1= mg/L
nan_inds=find(isnan(wq_conc_tmp));
not_nan_inds=setdiff(1:length(wq_conc_tmp),nan_inds)';
count=0;
for i=1:length(nan_inds)
    
    
    nan_ind_tmp=nan_inds(i);
    
    % find the nearest 30 samples that are not nans
    diff=[abs(not_nan_inds-nan_ind_tmp),not_nan_inds];
    sorted_diff=sortrows(diff);
    nearest_thirty=sort(sorted_diff(1:30,2));
    
    % read streamflow and wq data corresponding to these nearest 30 samples
    constituent.calib.C=wq_conc_tmp(nearest_thirty);
    constituent.calib.Q=strm(nearest_thirty);
    constituent.calib.date=date_wq_split(nearest_thirty,:);
    
    constituent.est.Q=strm(nan_ind_tmp);
    constituent.est.date=date_wq_split(nan_ind_tmp,:);
    
    % carry out bayesian regression
    res=bayesian_regress(constituent,modelNo,unitno);
    conc=res.est.unbiasedLoad/constituent.est.Q*(10^3/3600/24*35.3);        % estimated concentration in mg/L
    
    % fill in the missing data
    wq_conc_tmp(nan_ind_tmp)=conc;
    
end
wq_load=wq_conc_tmp.*strm*(10^(-3)/35.3);
% save reconstructed wq concentration data
sname='maumee_reconstructed_TP_conc.txt';
save_filename=fullfile(direc_htlp,sname);
fid=fopen(save_filename,'w');
fprintf(fid,'%s\t%s\t%s\n','Date','TP(mg/L)','TP_load(Kg/day)');
for write_ind=1:length(date)
    
    fprintf(fid,'%s\t%d\t%d\n',date{write_ind},wq_conc_tmp(write_ind),wq_load(write_ind));
    
end
fclose(fid);