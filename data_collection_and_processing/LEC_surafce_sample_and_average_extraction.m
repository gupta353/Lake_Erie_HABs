% This routine extracts surface data from LEC, and computes average of property
% of entire profile

clear all
close all
clc


direc=['D:/Research/EPA_Project/Lake_Erie_HAB/Data/lake_erie_LEC'];

station_list={'PhyAtt_4P.txt',...
    'PhyAtt_7M',...
    'PhyAtt_8M',...
    'PhyAtt_BUOY'...
    'PhyAtt_CRIB'...
    'PhyAtt_GR1',...
    'PhyAtt_MB18',...
    'PhyAtt_MB20',...
    'PhyAtt_OREGON'};

datenum_wrapper=@(x)datenum(x,'mm/dd/yyyy');
for station_ind=1:length(station_list)      % for loop for each station
    
    % read station data
    temp_station=station_list{station_ind};
    filename=fullfile(direc,temp_station);
    fid=fopen(filename,'r');
    formatspec=[repmat('%s',1,2),repmat('%f',1,35)];
    tline=fgetl(fid);
    data=textscan(fid,formatspec,'delimiter','\t');
    fclose(fid);
    varNames=strsplit(tline,'\t');
    
    date=data{1};
    time=data{2};
    depth=data{3};
    prop_data=cat(2,data{4:end});
    
    unique_dates=unique(date);
    
    for date_ind=1:length(unique_dates)     % for loop for date
        
        % find date, depth, and properties corresponding to this date
        ind=find(strcmp(date,unique_dates{date_ind}));
        temp_dates=date(ind);
        temp_times=time(ind);
        temp_depths=depth(ind);
        temp_prop_data=prop_data(ind,:);
        min_depth_ind=find(temp_depths==min(temp_depths));
        
        % prperties correspondig to surface depth and mean of properties at
        % all depths
        temp_mean_write_data=mean(temp_prop_data);
        temp_surface_write_data=temp_prop_data(min_depth_ind,:);
        surface_depth=temp_depths(min_depth_ind);
        mean_surface_depth=mean(temp_depths);
        
        % write surafce depth and mean data to cell arrays
        write_data_surface(date_ind,:)={temp_dates{1},temp_times{1},...
            [surface_depth,temp_surface_write_data]};
        write_data_mean(date_ind,:)={temp_dates{1},temp_times{1},...
            [mean_surface_depth,temp_mean_write_data]};
        
        
    end
    
    % sort the data according to date, for both surface and mean data
    date_surface=write_data_surface(:,1);
    datenum_surface=cellfun(datenum_wrapper,date_surface);
    datenum_surface(:,2)=(1:length(datenum_surface))';
    datenum_surface=sortrows(datenum_surface);
    write_data_surface=write_data_surface(datenum_surface(:,2),:);
    
    
    date_mean=write_data_mean(:,1);
    datenum_mean=cellfun(datenum_wrapper,date_mean);
    datenum_mean(:,2)=(1:length(datenum_mean))';
    datenum_mean=sortrows(datenum_mean);
    write_data_mean=write_data_mean(datenum_mean(:,2),:);
    
    % write data to text-files
    surface_fname=strcat('surface_',temp_station);
    surface_filename=fullfile(direc,surface_fname);
    wfid1=fopen(surface_filename,'wt');
    formatspec=[repmat('%s\t',1,36),'%s\n'];
    fprintf(wfid1,formatspec,varNames{:});
    
    mean_fname=strcat('mean_',temp_station);
    mean_filename=fullfile(direc,mean_fname);
    wfid2=fopen(mean_filename,'wt');
    formatspec=[repmat('%s\t',1,36),'%s\n'];
    fprintf(wfid2,formatspec,varNames{:});
    
    formatspec=['%s\t%s\t',repmat('%f\t',1,34),'%f\n'];
    for write_ind=1:length(write_data_surface)      % for loop to write data in textfiles, row-wise
        
        fprintf(wfid1,formatspec,write_data_surface{write_ind,1},...
            write_data_surface{write_ind,2},write_data_surface{write_ind,3});
        fprintf(wfid2,formatspec,write_data_mean{write_ind,1},...
            write_data_mean{write_ind,2},write_data_mean{write_ind,3});
    
    end
    fclose(wfid1); fclose(wfid2);
    
end