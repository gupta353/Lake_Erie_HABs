% This routine extracts surface data from NOAA, and computes average of
% property of profile

clear all
close all
clc

direc=['D:/Research/EPA_Project/Lake_Erie_HAB/Data/lake_erie_NOAA/'...
    'hab_field_sampling_and_nutrient_buoys/1.1/data/0-data'];

station_list={'hab_WE2.txt',...
    'hab_WE4.txt',...
    'hab_WE6.txt',...
    'hab_WE8.txt',...
    'hab_WE9.txt',...
    'hab_WE12.txt',...
    'hab_WE13.txt',...
    'hab_WE14.txt',...
    'hab_WE15.txt',...
    'hab_WE16.txt'};

datenum_wrapper=@(x)datenum(x,'mm/dd/yyyy');        % check the date format
for station_ind=1:length(station_list)      % for loop for each station
    
    % read data from textfile
    temp_station=station_list{station_ind};
    filename=fullfile(direc,temp_station);
    formatspec=repmat('%s',1,35);
    fid=fopen(filename,'r');
    tline=fgetl(fid);
    data=textscan(fid,formatspec,'delimiter','\t');
    fclose(fid);
    varNames=strsplit(tline,'\t');
    
    date=data{1};
    time=data{5};
    depth=data{3};
    depth_label=data{4};
    mid_data=cat(2,data{6:9});
    sky=data{10};                    % sky condition
    prop_data=cat(2,data{11:end});   % property data
    
    unique_dates=unique(date);
    
    for date_ind=1:length(unique_dates) % 'for' loop for each date
        
        % find data corresponding to a particular date
        ind=find(strcmp(date,unique_dates{date_ind}));
        temp_dates=date(ind);
        temp_times=time(ind);
        temp_depth=depth(ind);
        temp_depth_label=depth_label(ind);
        temp_mid_data=mid_data(ind,:);
        temp_sky=sky(ind);
        temp_prop_data=prop_data(ind,:);
        temp_depth_label_ind=find(strcmp(temp_depth_label,'Surface'));
        
        if isempty(temp_times(temp_depth_label_ind))        % in-case time-stamp is missing
            temp_times(temp_depth_label_ind)='missing';
        end
            
        write_data(date_ind,:)=[temp_dates{1},temp_times(temp_depth_label_ind),...
            temp_depth(temp_depth_label_ind),temp_sky(temp_depth_label_ind),...
            temp_mid_data(temp_depth_label_ind,:),temp_prop_data(temp_depth_label_ind,:)];
        
    end
    
    % sort data according to dates
    write_dates=write_data(:,1);
    write_datenums=cellfun(datenum_wrapper,write_dates);
    write_datenums(:,2)=1:length(write_datenums);
    write_datenums=sortrows(write_datenums);
    write_data=write_data(write_datenums(:,2),:);
    
    % write data to a text-file
    write_fname=strcat('surface_',temp_station);
    write_filename=fullfile(direc,write_fname);
    wfid=fopen(write_filename,'wt');
    formatspec=[repmat('%s\t',1,32),'%s\n'];
    fprintf(wfid,formatspec,varNames{1},varNames{5},varNames{3},...
        varNames{10},varNames{6:9},varNames{11:end});
    
    for write_ind=1:size(write_data,1)
        fprintf(wfid,formatspec,write_data{write_ind,:});
    end
    fclose(wfid);

end