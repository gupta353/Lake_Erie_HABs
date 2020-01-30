% routine to download streamflow data from USGS website

clear all
close all
clc
fclose('all');
%% path specification to download the data
path=['D:/Research/EPA_Project/Lake_Erie_HAB/Data'...       % path to download directory
    '/streamflow_data'];
%% load text file containing site numbers
%
fileID=fopen(['']);
site=textscan(fileID,'%s %s');
site=site{2};
n=size(site,1);                                     % number of sites in the file
%}
%% download the html pages
%
for i=1:n
    url_dummy=['https://waterdata.usgs.gov'...      % Sample url for relevant usgs website
        '/nwis/inventory/?site_no=03347000&agency_cd=USGS&amp'];
%   str_n=strcat('site_no=','0',num2str(site(i)));
    str_n=strcat('site_no=',site{i});
    url=regexprep(url_dummy,'site_no=03347000',str_n);
%     filename=num2str(site(i));
    filename=strcat(path,'/',site{i});
    openfile=urlwrite(url,filename);
end
%}
%% Read latitude, longitude, period of streamflow from HTML pages and streamflow data collection
%
mis_count=0;                                                % count for missing streamflow-data sites
avail_count=0;                                              % count for available streamflow-data sites
for j=1:n
    filename=fullfile(path,site{j});                        % file to be read    
    opened_file=fileread(filename);                         % read the file as character vector
    exp_loc=['Latitude\s*(?<lat>[^ ,]*),'...                % expression to match for lat-long
    '\s*\S+\s*?Longitude\s*(?<long>\S*)\s*'...
    '\S+\s*(?<proj>\w*)'];
    location{j}=regexp(opened_file,...                      % Cell-array of structures contaning location data
    exp_loc,'names');
    exp_period=['<td>&nbsp;(?<begin_date>\d+-\d+-\d+)&nbsp;'...
        '</td><td>&nbsp;(?<end_date>\d+-\d+-\d+)&nbsp;</td><td>&nbsp;</td>'];

    period{j}=regexp(opened_file,exp_period,'names');       % Cell-array of structures contaning record-length data
    % streamflow data download
    %
    if isempty(period{j})
        mis_count=mis_count+1;
        missing_sites(mis_count)=site(j);
    else
       url_dummy=['https://nwis.waterdata.usgs.gov/'...
           'il/nwis/uv/?cb_00060=on&cb_00065=on&'...
           'format=rdb&site_no=03335500&period=&begin'...
           '_date=1988-06-01&end_date=2006-12-31'];
       str_n=strcat('site_no=',site{j});                    % current site-number
       url=regexprep(url_dummy,'site_no=03335500',...       % replace dummy site number with current site number
            str_n);
       url=regexprep(url,'begin_date=1988-06-01',...        % replace dummy begin-date with current begin-date
            strcat('begin_date=',period{j}(1).begin_date));
       url=regexprep(url,'end_date=2006-12-31',...          % replace dummy end-date with current end-date
            strcat('end_date=',period{j}(1).end_date));
       loc=fullfile(path,strcat('streamflow_',site{j}));
       output=urlwrite(url,loc);                            % downloading streanflow data as a text file
       avail_count=avail_count+1;
       avail_sites(avail_count)=site(j);
    end
    %}
end
%}