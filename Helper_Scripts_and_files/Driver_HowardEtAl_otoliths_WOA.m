% Driver_HowardEtAl_otoliths_WOA.m
% This script and its associated functions forms part of the code
% for the analysis published in 
% Howard, E., at al.: Hypoxia traits imprinted in otolith delta13C
% from individual to global scales; Scientific Reports, 2024
%
% Author: H. Frenzel, UW-CICOES / NOAA-PMEL
%
% Date: August 27, 2024

data_dir = './WOA'; % on the local drive; adapt as needed
vnames_in = {'temperature'; 'salinity'; 'oxygen'; 'phosphate'};
WOA.vnames = {'temp','salt','o2','po4'};
fn_woa_combined = sprintf('%s/woa18_combined.nc', data_dir); % final output file

if ~exist(data_dir, 'dir')
    mkdir('.', data_dir)
end

% download one file from World Ocean Atlas 2009,
% used only for its vertical levels
woa09_base = 'https://www.ncei.noaa.gov/data/oceans/woa/WOA09/NetCDFdata/';
woa09_file = 'temperature_annual_1deg.nc';
local_file = sprintf('%s/%s', data_dir, woa09_file);
if ~exist(local_file, 'file')
    url = [woa09_base, woa09_file];
    websave(local_file, url);
end
woa09_depths = ncread(local_file, 'depth');

% download World Ocean Atlas 2018 data
nvars = length(WOA.vnames);
fn_woa_ann = cell(nvars,1);
short_var_name = cell(nvars,1);
period_name = cell(nvars,1);
fn_woa_combined_var = cell(nvars,1);
for v = 1:length(vnames_in)
    [fn_woa_ann{v}, short_var_name{v}, period_name{v}] = ...
        download_woa(vnames_in{v}, data_dir);
    fn_woa_combined_var{v} = combine_woa_data(short_var_name{v}, ...
        [data_dir, filesep, fn_woa_ann{v}]);
end

% create the starting point of the combined file
% ncks is far more efficient to use than native Matlab functions for this 
% task
cmd = sprintf('ncks -O -v lon,lat,time %s %s', fn_woa_combined_var{1}, ...
    fn_woa_combined);
system(cmd);

woa18_depths = ncread(fn_woa_combined_var{1}, 'depth');

% add the actual WOA fields to the combined file
create_woa_netcdf(fn_woa_combined, fn_woa_combined_var, ...
    short_var_name, WOA.vnames, woa09_depths, woa18_depths)

% read data from combined netcdf file
lon = ncread(fn_woa_combined, 'lon');
WOA.x = circshift(lon,180);
WOA.x(WOA.x < 0) = WOA.x(WOA.x < 0); % use lon = 0..360 degrees
WOA.y = ncread(fn_woa_combined, 'lat');
[~,~,area] = grid_area(WOA.x, WOA.y);
WOA.z = woa09_depths;
WOA.o2 = circshift(ncread(fn_woa_combined,'o2'), 180);
[nx,ny,nz,nt] = size(WOA.o2);
WOA.temp = circshift(ncread(fn_woa_combined,'temp'), 180);
WOA.salt = circshift(ncread(fn_woa_combined, 'salt'), 180);
WOA.po4 = circshift(ncread(fn_woa_combined, 'po4'), 180);
tmp = repmat(woa09_depths,[1,nx,ny,nt]);
z4d = permute(tmp,[2,3,1,4]);
WOA.po2 = O2pressure(WOA.o2,WOA.temp,WOA.salt,z4d);

% grid-related variables
landmask = WOA.o2(:,:,1,1) ./ WOA.o2(:,:,1,1);
area = area .* landmask;
WOA.a3d = repmat(area, [1 1 nz]);
WOA.zbot = 0.5 * (woa09_depths(1:end-1) + woa09_depths(2:end))';
WOA.zbot(end+1) = woa09_depths(end);
WOA.ztop(1) = 0;
WOA.ztop(2:length(WOA.zbot)) = WOA.zbot(1:end-1);
WOA.zedge = unique([WOA.zbot, WOA.ztop]);
WOA.dz = WOA.zbot - WOA.ztop;
dz3d = permute(repmat(WOA.dz', [1 nx ny]), [2 3 1]);
WOA.v3d = WOA.a3d .* dz3d;

% output to mat file
save('WOA.mat','WOA')