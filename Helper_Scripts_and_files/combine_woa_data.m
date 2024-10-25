function file_path_combined = combine_woa_data(short_var_name, file_path_ann)
% COMBINE_WOA_DATA This function combines data from monthly/seasonal/annual
% World Ocean Atlas files for one variable and creates a combined file.
% For layers that have monthly data, monthly values are used.
% For layers that do not have monthly, but seasonal values, seasonal
% values are linearly interpolated onto the months.
% For layers that have only annual values, these values are copied
% for each month.
%
% Author: H. Frenzel, UW-CICOES / NOAA-PMEL
%
% Date: August 27, 2024

var_name = [short_var_name, '_an']; % objectively analyzed

file_path_jan = strrep(file_path_ann, '00_01.nc', '01_01.nc');
file_path_win = strrep(file_path_ann, '00_01.nc', '13_01.nc');
file_path_combined = strrep(file_path_ann, '00_01.nc', '.nc');

vars_to_copy = {'lat'; 'lat_bnds'; 'lon'; 'lon_bnds'; 'depth'; ...
    'depth_bnds'; 'time'; var_name};

% Open the annual NetCDF file in read-only mode
ncid_ann = netcdf.open(file_path_ann, 'NC_NOWRITE');

% Create the combined NetCDF file
ncid_com = netcdf.create(file_path_combined, 'CLOBBER');

% Copy dimensions
[n_dims, ~, n_global_atts, ~] = netcdf.inq(ncid_ann);

for i = 0:n_dims-1
    [dim_name, dim_length] = netcdf.inqDim(ncid_ann, i);
    if strcmp(dim_name, 'time')
        dim_length = 12; % expand from annual to monthly data
    end
    netcdf.defDim(ncid_com, dim_name, dim_length);
end

% Copy global attributes
nc_glob = netcdf.getConstant('NC_GLOBAL');
for i = 0:n_global_atts-1
    gatt_name = netcdf.inqAttName(ncid_ann, nc_glob, i);
    gatt_value = netcdf.getAtt(ncid_ann, nc_glob, gatt_name);
    if strcmp(gatt_name, 'title')
        gatt_value = strrep(gatt_value, 'Annual', 'combined');
    end
    netcdf.putAtt(ncid_com, nc_glob, gatt_name, gatt_value);
end
% modify and create new global attributes
netcdf.putAtt(ncid_com, nc_glob, 'date_created', datestr(today))
netcdf.putAtt(ncid_com, nc_glob, 'date_modified', datestr(today))
netcdf.putAtt(ncid_com, nc_glob, 'history', ...
    'Created with combine_woa_data.m from individual files')

% Copy definitions of selected variables
varid_com = nan(length(vars_to_copy),1);
for v = 1:length(vars_to_copy)
    varid = netcdf.inqVarID(ncid_ann,vars_to_copy{v});
    [this_var_name, xtype, dim_IDs, natts] = netcdf.inqVar(ncid_ann, varid);
    varid_com(v) = netcdf.defVar(ncid_com, this_var_name, xtype, dim_IDs);
    if strcmp(this_var_name, 'time')
        varid_time_com = varid_com(v);
    end
    % Copy variable attributes
    for j = 0:natts-1
        att_name = netcdf.inqAttName(ncid_ann, varid, j);
        att_value = netcdf.getAtt(ncid_ann, varid, att_name);
        netcdf.putAtt(ncid_com, varid_com(v), att_name, att_value);
    end
end
netcdf.putAtt(ncid_com, varid_time_com, 'units', 'months')
netcdf.delAtt(ncid_com, varid_time_com, 'climatology')

netcdf.endDef(ncid_com)

% Copy values of selected variables (not time or the main variable)
for v = 1:length(vars_to_copy) - 2
    varid = netcdf.inqVarID(ncid_ann,vars_to_copy{v});
    %[this_var_name, xtype, dim_IDs, natts] = netcdf.inqVar(ncid_ann, varid);
    this_var = netcdf.getVar(ncid_ann, varid);
    netcdf.putVar(ncid_com, varid_com(v), this_var);
end
% time values are mid-month; change units to generic months
netcdf.putVar(ncid_com, varid_time_com, 0.5:11.5);

% determine the number of depth levels that are available in the
% different types of files (the annual file always has all depth levels)
depth_ann = ncread(file_path_ann, 'depth');
depth_mon = ncread(file_path_jan, 'depth');
depth_sea = ncread(file_path_win, 'depth');

nz_ann = length(depth_ann);
nz_mon = length(depth_mon);
nz_sea = length(depth_sea);

% type 1: T,S,O - nz_ann == nz_sea: 102; nz_mon: 57
% type 2: P     - nz_ann: 102; nz_sea == nz_mon: 43
if nz_ann == nz_sea && nz_sea > nz_mon
    type = 1;
elseif nz_ann > nz_sea && nz_sea == nz_mon
    type = 2;
else
    error('case not coded yet')
end

if type == 1
    % the deeper layers have seasonal data
    var_deep = [];
    for s = 1:4
        str_sea = sprintf('%02d_01', s+12);
        file_path_sea = strrep(file_path_jan, '01_01', str_sea);
        tmp = ncread(file_path_sea, var_name, [1 1 nz_mon+1 1], ...
            [inf inf nz_sea-nz_mon 1]); % 1-offset
        [nx_sea, ny_sea, ~] = size(tmp);
        var_deep = cat(4,var_deep,tmp);
    end
    for mth = 1:12
        if mod(mth, 3) > 0
            season = 1 + floor(mth / 3); % "main" season
        else
            season = mth / 3;
        end
        if mod(mth, 3) == 2
            % month is in middle of season: copy seasonal values
            netcdf.putVar(ncid_com, varid_com(end), [0 0 nz_mon mth-1], ...
                [nx_sea ny_sea nz_sea-nz_mon 1], var_deep(:,:,:,season)); % 0-offset
        else
            if mod(mth, 3) == 1 % one month before middle of season
                season2 = season - 1;
                if season2 == 0
                    season2 = 4;
                end
            else % one month after middle of season
                season2 = season + 1;
                if season2 == 5
                    season2 = 1;
                end
            end
            this_var_deep = (2.0 * var_deep(:,:,:,season) + ...
                var_deep(:,:,:,season2)) / 3.0;
            netcdf.putVar(ncid_com, varid_com(end), [0 0 nz_mon mth-1], ...
                [nx_sea ny_sea nz_sea-nz_mon 1], this_var_deep); % 0-offset
        end
    end
elseif type == 2
    % the deeper layers have only annual data
    var_deep = ncread(file_path_ann, var_name, [1 1 nz_mon+1 1], ...
        [inf inf nz_ann-nz_mon 1]); % 1-offset
    [nx_ann, ny_ann, ~] = size(var_deep);
    for mth = 1:12
        netcdf.putVar(ncid_com, varid_com(end), [0 0 nz_mon mth-1], ...
            [nx_ann ny_ann nz_ann-nz_mon 1], var_deep); % 0-offset
    end
end
nx = size(var_deep, 1);
ny = size(var_deep, 2);

for mth = 1:12
    % shallow layers: copy monthly data
    str_mth = sprintf('%02d_01', mth);
    file_path_mon = strrep(file_path_jan, '01_01', str_mth);
    var_shallow = ncread(file_path_mon, var_name, [1 1 1 1], ...
        [inf inf nz_mon 1]); % 1-offset
    netcdf.putVar(ncid_com, varid_com(end), [0 0 0 mth-1], ...
        [nx ny nz_mon 1], var_shallow); % 0-offset
end

% Close both NetCDF files
netcdf.close(ncid_ann);
netcdf.close(ncid_com);
