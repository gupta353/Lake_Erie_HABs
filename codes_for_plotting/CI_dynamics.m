

clear all
close all
clc

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/';
agency = 'lake_erie_NOAA';
station = 'erie-cmt_bydate';
met_par = 'daily_wspeed.txt';

% read CI data
fname = 'CI_total_combined_MERIS.txt';
filename = fullfile(direc,'remote_sensing_data',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);

CI_dates = data{1};
CI_vals = data{2};
CI_datenums = cellfun(@(x)datenum(x,'yyyy-mm-dd'),CI_dates);

%%
years = {'2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012'};
years_datenum = cellfun(@(x)datenum(x,'yyyy'),years);
months = {'05','06','07','08','09','10'};


for year_ind = 1:length(years_datenum)-1
    
    ind1 = find(CI_datenums>=years_datenum(year_ind)& CI_datenums<years_datenum(year_ind+1));
    comp_datenums = (CI_datenums(ind1(1)):10:CI_datenums(ind1(end)))';
    CI_tmp = NaN*ones(length(comp_datenums),1);
    
    for ii = 1:length(ind1)
        ind = ind1(ii);
        CI_tmp(comp_datenums==CI_datenums(ind)) = CI_vals(ind);
    end
    
    max_CI = max(CI_tmp);
    plot(CI_tmp(1:end-1)/max_CI,CI_tmp(2:end)/max_CI,'linewidth',2); hold on
    plot(CI_vals(ind1(1:end-1))/max_CI,CI_vals(ind1(2:end))/max_CI,'--')
    scatter(CI_tmp(1)/max_CI,CI_tmp(2)/max_CI,'facecolor','b')
    scatter(CI_tmp(end-1)/max_CI,CI_tmp(end)/max_CI,'facecolor','r')
    title(['year = ',years{year_ind},', Maximum CI = ',num2str(max_CI)],'fontname','arial','fontsize',12)
    xlabel('Relative CI(t)','fontname','arial','fontsize',12)
    ylabel('Relative CI(t+1)','fontname','arial','fontsize',12);
    set(gca,'fontname','arial','fontsize',12);
    
    sname = ['dynamics_of__CI_',years{year_ind},'.svg'];
    filename = fullfile('D:/Research/EPA_Project/Lake_Erie_HAB','matlab_codes/plots',sname);
    %         print(filename,'-r300','-djpeg');
    fig2svg(filename);
    pause;
    close all
end





