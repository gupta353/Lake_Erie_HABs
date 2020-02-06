% this routine processes data from the Lake Erie NOAA center, and arranges
% data corresponding to each station in a single text-files with
% meterological data

clear all
close all
clc

direc=['D:/Research/EPA_Project/Lake_Erie_HAB/Data/lake_erie_NOAA'];
format_dir='station_data_format';           % name of the folder which contains station format data

% format files for different station files
station_files={'erie-clv_bydate+data'};

% names of the variables that are needed to be read and arranged
variable_names={'atemp(DegC)'};

% what kind of data is to be extracted - met, temp or ysi
data_tobeextracted='data';

for station_ind=1%:length(station_files)                 % 'for' loop for each station
    
    % read station format files
    file_contain_format_data=station_files{station_ind};
    file_contain_format_data_split=strsplit(file_contain_format_data,'+');
    folder_file_name=file_contain_format_data_split{1};        % name of folder (whcih contains files)
    filename=fullfile(direc,format_dir,...
        strcat(file_contain_format_data,'.txt'));
    fid=fopen(filename,'r');
    data=textscan(fid,'%s%s%d%s%s','delimiter','\t','headerlines',1);
    fclose(fid);
    
    folder_name=data{1};
    year=data{2};
    num_columns_toberead=data{3};
    cols_toberead=data{4};
    colnames_toberead=data{5};
    
    % create text files in which data will be written
    write_filename_atemp=fullfile(direc,folder_file_name,'atemp_2003.txt');

%     write_filename_wtemp=fullfile(direc,folder_file_name,'wtemp.txt');
    fid_atemp=fopen(write_filename_atemp,'w');
    
    fprintf(fid_atemp,[repmat('%s\t',1,2),'%s\n'],'Date(MM/dd/yyyy)','time(hh:mm:ss)',variable_names{1:end});
    
    % write data corresponding the station 'folder_file_name' in four text file
    for year_ind=1:length(year)                   % 'for' loop for each year in a station
        
        temp_year=year{year_ind};
        temp_folder_name=folder_name{year_ind};
        colsequence_toberead_split=strsplit(cols_toberead{year_ind},',');
        colsequence_toberead_split=cellfun(@str2double,colsequence_toberead_split);
        colnames_toberaed_split=strsplit(colnames_toberead{year_ind},',');
        temp_num_columns_toberead=num_columns_toberead(year_ind);
        
        % list all the files in the folder containing textfiles
        list=dir(fullfile(direc,folder_file_name,temp_year));
        list_names={list.name};
        
        if sum(strcmp(list_names,'data'))==1                 % if 'data' files are contained in the folder
            
            data_type_folder=fullfile(direc,folder_file_name,temp_year,data_tobeextracted); % folder containing three data types - met, temp, ysi
            list_textfiles=dir(data_type_folder);
            
            for textfile_ind=1:size(list_textfiles,1)       % 'for' loop for each textfile in a year
                
                textfile_name=list_textfiles(textfile_ind).name;
                
                if ~strcmp(textfile_name,{'.','..'})
                    
                    fid_textfile=fopen(fullfile(data_type_folder,textfile_name),'r');
                    data=textscan(fid_textfile,repmat('%s',1,temp_num_columns_toberead),...
                        'delimiter',',','headerlines',1);
                    fclose(fid_textfile);
                    
                    
                    % read atemp, wspeed, wdir, and wtemp data in separate
                    % variables
                    data_tobewritten=data(colsequence_toberead_split);
                    data_tobewritten_atemp=cat(2,data_tobewritten{:});
                    data_tobewritten_atemp_date_time=data_tobewritten_atemp(:,1:2);
                    data_tobewritten_atemp_param=cellfun(@str2double,data_tobewritten_atemp(:,3));
                    
                    % write data in a 'for' loop
                    for writedatarow_ind=1:size(data_tobewritten_atemp,1)
                        fprintf(fid_atemp,'%s\t%s\t%f\n',...
                            data_tobewritten_atemp_date_time{writedatarow_ind,1:2},...
                            data_tobewritten_atemp_param(writedatarow_ind,:));
                    end
                    
                end
                
            end
            
        end
    end
    fclose(fid_atemp);
    
end
