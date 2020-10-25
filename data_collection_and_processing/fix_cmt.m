% fix format 'erie-cmt_by_date' data

clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/lake_erie_NOAA/meteorological_data/erie-cmt_bydate';
year='2009';
mkdir(fullfile(direc,year,'metcsi_updated'));

list = dir(fullfile(direc,year,'metcsi_raw'));
for fname_ind = 1:size(list,1)
    fname = list(fname_ind).name;
    if ~strcmp(fname,{'.','..'})
        filename = fullfile(direc,year,'metcsi_raw',fname);
        fid = fopen(filename,'r');
        formatspec = [repmat('%s',1,1),'%s'];
        headerline = fgetl(fid);
        data = textscan(fid,formatspec);
        fclose(fid);
        
        formatted_data = {};
        for dind = 1:length(data{1})
            
            tmp_data = data{1}{dind};
            split_data = strsplit(tmp_data,'""');
            
            if length(split_data)<=4
                formatted_data = [formatted_data;tmp_data];
            else
                formatted_data = [formatted_data;strcat(split_data{1},'""',split_data{2},'","",,,,,,,,,,,,');strcat('"',split_data{3},'""',split_data{4})];
            end
            
        end
        
        % write data to a textfile with same name as of original file but in a
        % different folder
        
        filename = fullfile(direc,year,'metcsi_updated',fname);
        fid = fopen(filename,'w');
        fprintf(fid,'%s\n',headerline);
        for wind = 1:length(formatted_data)
            fprintf(fid,'%s\n',formatted_data{wind});
        end
        fclose(fid);
    end
end