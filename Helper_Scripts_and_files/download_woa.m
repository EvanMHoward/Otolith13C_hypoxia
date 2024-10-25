function [file_name_ann, short_name, period_name] = ...
    download_woa(var_name, data_dir)
% DOWNLOAD_WOA This function downloads all monthly, seasonal, and annual
% files for one specified variables from the World Ocean Atlas.
%
% Author: H. Frenzel, UW-CICOES / NOAA-PMEL
%
% Date: August 27, 2024


base_url = 'https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/';

if strcmp(var_name , 'temperature')
    short_name = 't';
    period_name = 'decav';
elseif strcmp(var_name , 'salinity')
    short_name = 's';
    period_name = 'decav';
elseif strcmp(var_name , 'oxygen')
    short_name = 'o';
    period_name = 'all';
elseif strcmp(var_name , 'phosphate')
    short_name = 'p';
    period_name = 'all';
else
    error('variable not yet coded')
end

% 0 is annual, 1 to 12 monthly, 13 to 16 seasonal
% monthly values are not available for the deep ocean
for time_period = 0:16
    file_name = sprintf('woa18_%s_%s%02d_01.nc', period_name, ...
        short_name, time_period);
    if time_period == 0
        file_name_ann = file_name;
    end
    url = sprintf('%s%s/netcdf/%s/1.00/%s', base_url, var_name, ...
        period_name, file_name);
    local_file = sprintf('%s/%s', data_dir, file_name);
    if ~exist(local_file, 'file')
        websave(local_file, url);
    end
end
