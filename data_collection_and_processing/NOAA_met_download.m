% this routime downloads meteorological data over Lake Erie from NOAA
% archive (https://www.glerl.noaa.gov/res/recon/data/) This code work in
% later versions of MATLAB (not in 2014a) 
% README files have the following exceptions which needs to be downloaded
% separately: 
%       0-REDME-2010-clv-met-note.txt;
%       0-REDME-2010-clv-met-ERRORS.txt;
%       0-REDME-2010-cln-met.txt
%       https://www.glerl.noaa.gov/res/recon/data/bystation/erie-clv/2013/ysi/ysi20130503.dat
%       all erie_cm_bydate-2004 data has to be downloaded manually
clear all
close all
clc

direc=['D:/Research/EPA_Project/Lake_Erie_HAB/Data/lake_erie_NOAA/meteorological_data'];

metadata_files={'download_metadata_erie-clv_bydate',...
    'download_metadata_erie-clv_byparam',...
    'download_metadata_erie-clvnorth_bydate',...
    'download_metadata_erie-cmt_bydate',...
    'download_metadata_erie-cmt_csv',...
    'download_metadata_erie-rbs_csv',...
    'download_metadata_erie-wle_bydate',...
    'download_metadata_erie-cmt_bydate_2016'};

partial_url='https://www.glerl.noaa.gov/res/recon/data/bystation';

% read each metadata file and download data
for metadata_fileno=8%:length(metadata_files)
    
    % read metafile
    name=metadata_files{metadata_fileno};
    str_split=strsplit(name,'_');
    station=str_split{3}; 
    archive_file_format=str_split{4};
    
    cd(direc)
    mkdir(strcat(station,'_',archive_file_format));
    cd(strcat(station,'_',archive_file_format))
    
    % download data
    if strcmp(archive_file_format,'bydate')                 % if data files are arranged date-wise
        
        filename=fullfile(direc,'meteorological_data_and_webcams_server_files_format',strcat(name,'.txt'));
        fid=fopen(filename,'r');
        data=textscan(fid,'%s%s%s%s%s%s%d','delimiter','\t');
        fclose(fid);
        
        year=data{1}; folder=data{2};
        begin_date=data{3}; end_date=data{4};
        filename_format=data{5}; filename_extension=data{6};
        fileformat_col_num=data{7};
        wrapper=@(x)datestr(x,'yyyymmdd');
        
        
        for year_fileno=1:length(year)
            
            temp_year=year{year_fileno};
            mkdir(temp_year)
            
            % read variabls for the particular row in 'download-metadata' file
            name_format=filename_format{year_fileno};
            name_extension=filename_extension{year_fileno};
            begin_datenum=datenum(begin_date(year_fileno),'yyyymmdd');
            end_datenum=datenum(end_date(year_fileno),'yyyymmdd');
            datenum_list=begin_datenum:end_datenum;
            datelist=datestr(datenum_list,'yyyymmdd');
            col_num= fileformat_col_num(year_fileno);
            
            for datelist_num=1:size(datelist,1)
                
                download_filename=strrep(name_format,'[date]',datelist(datelist_num,:));
                download_filename=strcat(download_filename,'.',name_extension);
                
                
                url=strcat(partial_url,'/',station,'/',temp_year,'/',...
                    folder{year_fileno},'/',download_filename);
                
                % using webread function
                myreadtable = @(filename)readtable(filename,...
                    'Format',strcat(repmat('%s',1,col_num)),...
                    'Delimiter','space','MultipleDelimsAsOne',1);
                options = weboptions('ContentReader',myreadtable);
                readdata=webread(url,options);
                write_name=strcat(temp_year,'/',download_filename);
                writetable(readdata,write_name);

                
            end
        end
        
    elseif strcmp(archive_file_format,'byparam')            % if data files are arranged parameter-wise
        
        filename=fullfile(direc,strcat(name,'.txt'));
        fid=fopen(filename,'r');
        data=textscan(fid,'%s%s%s%s%d','delimiter','\t');
        fclose(fid);
        
        year=data{1}; folder=data{2};
        param_list=data{3}; filename_extension=data{4};
        fileformat_col_num=data{5};
        
        for param_ind=1:length(param_list)
            
            param_name=param_list{param_ind};
            string_num=fileformat_col_num(param_ind);
            col_num= fileformat_col_num(param_ind);
            
            url=strcat(partial_url,'/',station,'/',year{param_ind},'/',...
                folder{param_ind},...
                '/',param_name,'.',filename_extension{param_ind});
            myreadtable = @(filename)readtable(filename,...
                'Format',strcat(repmat('%s',1,col_num)),...
                'Delimiter','space','MultipleDelimsAsOne',1);
            options = weboptions('ContentReader',myreadtable);
            readdata=webread(url,options);
            write_name=param_name;
            writetable(readdata,write_name);
        end
        
    elseif strcmp(archive_file_format,'csv')
        
        filename=fullfile(direc,strcat(name,'.txt'));
        fid=fopen(filename,'r');
        data=textscan(fid,'%s%s%s%s','delimiter','\t');
        fclose(fid);
        
        year=data{1}; folder=data{2};
        param_list=data{3}; filename_extension=data{4};
        
        for param_ind=1:length(param_list)
            
            param_name=param_list{param_ind};
            
            url=strcat(partial_url,'/',station,'/',year{param_ind},'/',...
                folder{param_ind},'/',param_name,'.',filename_extension{param_ind});
            write_name=strcat(param_name,'_',year{param_ind});
            options=weboptions('Timeout',100);
            readdata=websave(write_name,url,options);
        end
    end
    
    
end