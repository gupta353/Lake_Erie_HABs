% This routine computes flow-weighted-average-daily-scale concentrations at HTLP
% (Heidelberg Tributary loading program) stations
% Note: (1) all the samples with time-window greater than 1 are removed in this
% script

clear all
close all
clc

rows_tobe_removed=[793:796,815:818,7219:7221,8522:8524,9611:9613,...
    13558:13560,13961:13963,13970:13972,13979:13981,13988:13990,...
    13994:13996,14012:14014,14021:14023,14038:14040,14088:14090,...
    14106:14109,17393:17395,17398:17400,17603:17605];
direc='D:/Research/EPA_Project/Lake_Erie_HAB/Data/HTLP';

% variables to be extracted
var_names={'Datetime_dateAndTimeOfSampleCollection_',...
    'SampleTimeWindow_Days',...
    'Flow_CFS',...
    'SS_Mg_L_suspendedSolids_',...
    'TP_Mg_LAsP',...
    'SRP_Mg_L_AsP'...
    'NO23_Mg_LAsN',...
    'TKN_Mg_L_TotalKjeldahlNitrogen_'};

station_list={'maumeedata.txt'};

datesplit_wrapper=@(x)strsplit(x,char(32));
datenum_wrapper=@(x)datenum(x,'mm/dd/yyyy HH:MM');
for station_ind=1:length(station_list)          % read data contained in excel files one-by-one
    
    station_name=station_list{station_ind};
    
    if ~strcmp(station_name,'basic-station-information-ao-wy2019.xlsx')
        
        % read data
        filename=fullfile(direc,station_name);
        T=readtable(filename,'delimiter','\t');
        
        time_window=T.(var_names{2});
        ind_tw=find(time_window<=1);                                % the samples with time-window greater than 1
        ind_tw=setdiff(ind_tw,rows_tobe_removed);
        
        % datenums
        datetime=T.(var_names{1});
        
        % streamflow and other properties
        strm=T.(var_names{3});
        strm=cellfun(@str2double,strm);
        % find indices that correspond to negative streamflow
        ind_strm=find(strm>0);
        ind_tw=intersect(ind_tw,ind_strm);
        strm=strm(ind_tw)*0.028316847;   % streamflow in cms
        % selec sane datenums
         datetime=datetime(ind_tw);
        datenums=cellfun(datenum_wrapper,datetime);
        datenum_days=floor(datenums);
        % select sane time-window data
        time_window=time_window(ind_tw);
        unique_datenum_days=unique(datenum_days);
        % read property data
        SS=T.(var_names{4});        SS=SS(ind_tw);
        TP=T.(var_names{5});        TP=TP(ind_tw);
        SRP=T.(var_names{6});       SRP=SRP(ind_tw);
        NO3_NO2=T.(var_names{7});   NO3_NO2=NO3_NO2(ind_tw);
        TKN=T.(var_names{8});       TKN=TKN(ind_tw);
        
        conc=[SS,TP,SRP,NO3_NO2,TKN];
        
        % compute loads
        loads=bsxfun(@times,conc,strm)*10^(-3); % loads in Kg/s
        
         next_day_initial_time_window=[];
         next_day_initial_load=[];
         next_day_initial_strm=[];
         old_temp_date=0;
         daily_conc=NaN(length(unique_datenum_days),5);
         daily_flow=NaN(length(unique_datenum_days),1);
        % process data corresponding to each unique date in excel file in 'for' loop
        for date_ind=1:length(unique_datenum_days)      
            
            temp_date=unique_datenum_days(date_ind);
            
            if temp_date-1~=old_temp_date
                
                next_day_initial_time_window=[];
                next_day_initial_load=[];
                next_day_initial_strm=[];
                
            end
            
            ind=find(datenum_days==temp_date);
            temp_time_window=[next_day_initial_time_window;time_window(ind)];
            temp_loads=[next_day_initial_load;loads(ind,:)];
            temp_strm=[next_day_initial_strm;strm(ind)];
            
            if sum(temp_time_window)>1       % if total time-window is greater than one day
                
                % move a portion of teh time-window (which is greater than
                % 1) to next day
                last_time_window=1-sum(temp_time_window(1:end-1));
                next_day_initial_time_window=sum(temp_time_window)-1;
                next_day_initial_load=temp_loads(end,:);
                next_day_initial_strm=temp_strm(end);
                temp_time_window(end)=last_time_window;
                
            else                             % if total time-window is not greater than one day
                next_day_initial_time_window=[];
                next_day_initial_load=[];
                next_day_initial_strm=[];
            end
            
            temp_total_load=bsxfun(@times,temp_loads,temp_time_window)*24*3600;       % total load in Kg in each time-window
            temp_total_flow=temp_strm.*temp_time_window*24*3600;                      % total streamflow in cms in each time-window
            daily_conc(date_ind,:)=sum(temp_total_load,1)/sum(temp_total_flow)*10^3;  % average-daily concentration
            daily_flow(date_ind)=sum(temp_total_flow)/sum(temp_time_window)/24/3600;  % mean-daily flow
            
            if sum(temp_time_window<0)>0
                error('negative time window detected');
            end
            
        end
    end
    
    % save daily concentration data to a textfile
    sname=strcat('daily_',station_name);
    save_filename=fullfile(direc,sname);
    fid=fopen(save_filename,'w');
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n',var_names{1},'Flow_CMS',var_names{4:end});
    for write_ind=1:size(daily_conc,1)
        fprintf(fid,'%s\t%f\t%f\t%f\t%f\t%f\t%f\n',...
            datestr(unique_datenum_days(write_ind),'mm/dd/yyyy'),...
            daily_flow(write_ind),daily_conc(write_ind,:));
    end
    fclose(fid);
end