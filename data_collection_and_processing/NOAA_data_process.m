% this routine processes data from the Lake Erie NOAA center, and arranges
% data corresponding to each station in a single text-files with
% meterological data

clear all
close all
clc

direc=['D:/Research/EPA_Project/Lake_Erie_HAB/Data/lake_erie_NOAA'];
format_dir='station_data_format';           % name of the folder which contains station format data

% format files for different station files
station_files={'erie-clv_bydate+met',...
    'erie-clvnorth_bydate+met',...
    'erie-cmt_bydate+met',...
    'erie-wle_bydate+met'};

% names of the variables that are needed to be read and arranged
variables_names={'atemp(DegC)',...
    'wspeed_avg(m/s)',...
    'wdir_avg(Deg)',...
    'wtemp(DegC)'};

% what kind of data is to be extracted - met, temp or ysi
data_tobeextracted='met';

for station_ind=4%:length(station_files)                 % 'for' loop for each station
    
    % read station format files
    file_contain_format_data=station_files{station_ind};
    file_contain_format_data_split=strsplit(file_contain_format_data,'+');
    folder_file_name=file_contain_format_data_split{1};        % name of file (which conatines format info) and folder (whcih contains files)
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
    write_filename_atemp=fullfile(direc,folder_file_name,'atemp.txt');
    write_filename_wspeed=fullfile(direc,folder_file_name,'wspped_avg.txt');
    write_filename_wdir=fullfile(direc,folder_file_name,'wdir_avg.txt');
    write_filename_wtemp=fullfile(direc,folder_file_name,'wtemp.txt');
    fid_atemp=fopen(write_filename_atemp,'w');
    fid_wspeed=fopen(write_filename_wspeed,'w');
    fid_wdir=fopen(write_filename_wdir,'w');
    fid_wtemp=fopen(write_filename_wtemp,'w');
    fprintf(fid_atemp,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
        'year','month','day','hour','minute','sec','atemp(DegC)');
    fprintf(fid_wspeed,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
        'year','month','day','hour','minute','sec','wspeed_avg(m/s)');
    fprintf(fid_wdir,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
        'year','month','day','hour','minute','sec','wdir_avg(Deg)');
    fprintf(fid_wtemp,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
        'year','month','day','hour','minute','sec','wtemp(DegC)');
    
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
        
        if sum(strcmp(list_names,'met'))==1                 % if 'met' files are contained in the folder
            
            data_type_folder=fullfile(direc,folder_file_name,temp_year,'met'); % folder containing three data types - met, temp, ysi
            list_textfiles=dir(data_type_folder);
            
            for textfile_ind=3:size(list_textfiles,1)       % 'for' loop for each textfile in a year
                
                textfile_name=list_textfiles(textfile_ind).name;
                if ~strcmp(textfile_name,{'.','..'})
                    
                    fid_textfile=fopen(fullfile(data_type_folder,textfile_name),'r');
                    data=textscan(fid_textfile,repmat('%s',1,temp_num_columns_toberead),...
                        'delimiter',',','headerlines',1);
                    fclose(fid_textfile);
                    
                    
                    % read atemp, wspeed, wdir, and wtemp data in separate
                    % variables
                    data_tobewritten=data(colsequence_toberead_split);
                    data_tobewritten_atemp=cat(2,data_tobewritten{[1:6,7]});
                    data_tobewritten_atemp=cellfun(@str2double,data_tobewritten_atemp);
                    data_tobewritten_wspeed=cat(2,data_tobewritten{[1:6,8]});
                    data_tobewritten_wspeed=cellfun(@str2double,data_tobewritten_wspeed);
                    data_tobewritten_wdir=cat(2,data_tobewritten{[1:6,9]});
                    data_tobewritten_wdir=cellfun(@str2double,data_tobewritten_wdir);
                    data_tobewritten_wtemp=cat(2,data_tobewritten{[1:6,10]});
                    data_tobewritten_wtemp=cellfun(@str2double,data_tobewritten_wtemp);
                    
                    fprintf(fid_atemp,'%d\t%d\t%d\t%d\t%d\t%d\t%d\n',...
                        data_tobewritten_atemp');
                    fprintf(fid_wspeed,'%d\t%d\t%d\t%d\t%d\t%d\t%d\n',...
                        data_tobewritten_wspeed');
                    fprintf(fid_wdir,'%d\t%d\t%d\t%d\t%d\t%d\t%d\n',...
                        data_tobewritten_wdir');
                    fprintf(fid_wtemp,'%d\t%d\t%d\t%d\t%d\t%d\t%d\n',...
                        data_tobewritten_wtemp');
                    
                end
                
            end
            
        end
    end
    fclose(fid_atemp); fclose(fid_wspeed); fclose(fid_wdir); fclose(fid_wtemp);
    
end