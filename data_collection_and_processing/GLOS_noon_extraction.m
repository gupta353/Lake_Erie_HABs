% This routine extracts data from GLOS such that for each day the data
% closest to noon is collected

clear all
close all
clc

direc=['D:/Research/EPA_Project/Lake_Erie_HAB/Data/'...
    'lake_erie_GLOS/extracted_data'];
list=dir(direc);

strsplit_wrapper=@(x)strsplit(x,'T');
time_wrapper=@(x)datenum(x,'HH:MM:SS');
noon_num=time_wrapper('12:00:00');
for station_ind=10:size(list,1)         % 'for' loop for each station
    
    temp_station=list(station_ind).name;
    if ~strcmp(temp_station,{'.','..'})
        
        % read data
        filename=fullfile(direc,temp_station);
        fid=fopen(filename,'r');
        tline1=fgetl(fid);
        tline2=fgetl(fid);
        formatspec=['%s',repmat('%f',1,14)];
        data=textscan(fid,formatspec,'delimiter','\t');
        fclose(fid);
        varNames=strsplit(tline1,'\t');
        varUnits=strsplit(tline2,'\t');
        
        date_time=data{1};
        prop_data=cat(2,data{2:end});
        
        for split_ind=1:length(date_time)
            date_time_split=strsplit(date_time{split_ind},'T');
            date{split_ind}=date_time_split{1};
            time{split_ind}=date_time_split{2}(1:end-1);        % removal of 'Z' at the end
        end
        
        unique_date=unique(date);
        for date_ind=1:length(unique_date)      % 'for' loop for each date
            
            ind=find(strcmp(date,unique_date(date_ind)));
            temp_dates=date(ind);
            temp_times=time(ind);
            temp_times_nums=cellfun(time_wrapper,temp_times);
            temp_prop_data=prop_data(ind,:);
            
            % find a time closes to noon
            diff_num=abs(temp_times_nums-noon_num);
            min_ind=find(diff_num==min(diff_num));
            min_ind=min_ind(1);     % in-case of two min_inds, the former one is selected
            
            write_data(date_ind,:)={temp_dates{1},temp_times{min_ind},...
                temp_prop_data(min_ind,:)};
        end
        
        % write data to a text-file
        write_fname=strcat('noon_',temp_station);
        write_filename=fullfile(direc,write_fname);
        wfid=fopen(write_filename,'wt');
        formatspec=[repmat('%s\t',1,15),'%s\n'];
        fprintf(wfid,formatspec,'date',varNames{1:end});
        fprintf(wfid,formatspec,'UTC',varUnits{1:end});
        
        formatspec=['%s\t%s\t',repmat('%f\t',1,13),'%f\n'];
        for write_ind=1:size(write_data,1)      % write data in 'for' loop
            fprintf(wfid,formatspec,write_data{write_ind,1},...
                write_data{write_ind,2},write_data{write_ind,3});
        end
        fclose(wfid);
    end
    unique_date=[]; data=[]; date=[]; time=[]; prop_data=[];
    
    
end