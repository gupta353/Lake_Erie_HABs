% This routine extracts physical and biological data from
% GLOS-local-repository

clear all
close all
clc

direc='D:/Research/EPA_Project/Lake_Erie_HAB/Data/lake_erie_GLOS/';

propt={'time',...
    'Air_Temperature',...
    'Water_Temperature_at_Surface',...
    'Wind_Speed',...
    'Wind_from_Direction',...
    'dissolved_oxygen',...
    'dissolved_oxygen_saturation',...
    'pH',...
    'water_conductivity',...
    'ysi_blue_green_algae',...
    'ysi_chlorophyll',...
    'ysi_turbidity',...
    'depth',...
    'latitude',...
    'longitude'};

fname='station_format.txt';
filename=fullfile(direc,'station_format',fname);
fid=fopen(filename,'r');
data=textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);
station_names=data{1};
cols_num=data{2};

for station_ind=1:length(station_names)
    
    temp_station_name=station_names{station_ind};
    temp_col_num=cols_num(2);
    
    % read csv file
    temp_filename=fullfile(direc,'raw_data',...
        strcat(temp_station_name,'.csv'));
    data=readtable(temp_filename);
    varnames=data.Properties.VariableNames;
    
    
    % create a texfile to write data
    write_fname=strcat(temp_station_name,'.txt');
    write_filename=fullfile(direc,'extracted_data',write_fname);
    fid=fopen(write_filename,'wt');
%     formatspec=repmat('%s\t',1,length(propt));
%     fprintf(fid,formatspec,propt{1:end});
    
    % write data to another table
    write_data=table;
    for varname_ind=1:length(propt)
        
        temp_propt=propt{varname_ind};
        ind=find(strcmp(varnames,temp_propt));
        if isempty(ind)
            temp_table=table(NaN*ones(size(data,1),1),...
                'VariableNames',{temp_propt});
            write_data=[write_data,temp_table];
        else
            write_data=[write_data,data(:,ind)];
        end
    end
    writetable(write_data,write_filename,'Delimiter','\t');
    
end

