function create_woa_netcdf(fn_woa_combined, fn_woa_combined_var, ...
    short_var_name, WOA_var_out, woa09_depths, woa18_depths)
% CREATE_WOA_NETCDF This function reads multiple netcdf files with
% different variables and creates a combined netcdf file that contains
% all these variables.
% 
% Author: H. Frenzel, UW-CICOES / NOAA-PMEL
%
% Date: August 27, 2024 

ncid_out = netcdf.open(fn_woa_combined, 'WRITE');

dimids(1) = netcdf.inqDimID(ncid_out, 'lon');
dimids(2) = netcdf.inqDimID(ncid_out, 'lat');
dimids(4) = netcdf.inqDimID(ncid_out, 'time');

netcdf.reDef(ncid_out)
% define depth variable
dimids(3) = netcdf.defDim(ncid_out, 'depth', length(woa09_depths));
varid_z = netcdf.defVar(ncid_out, 'depth', 'NC_FLOAT', dimids(3));
netcdf.putAtt(ncid_out, varid_z, 'long_name', 'depth');
netcdf.putAtt(ncid_out, varid_z, 'standard_name', 'depth');
netcdf.putAtt(ncid_out, varid_z, 'units', 'm');
netcdf.putAtt(ncid_out, varid_z, 'positive', 'down');
netcdf.putAtt(ncid_out, varid_z, 'long_name', 'depth');
netcdf.putAtt(ncid_out, varid_z, 'axis', 'Z');
netcdf.putAtt(ncid_out, varid_z, 'description', 'Standard Depth Levels');

varid_woa = cell(length(WOA_var_out));
var_name_in = cell(length(WOA_var_out));
for v = 1:length(WOA_var_out)
    varid_woa{v} = netcdf.defVar(ncid_out, WOA_var_out{v}, 'NC_FLOAT', dimids);
    ncid_in = netcdf.open(fn_woa_combined_var{v}, 'NOWRITE');
    var_name_in{v} = [short_var_name{v}, '_an'];
    varid_in = netcdf.inqVarID(ncid_in, var_name_in{v});
    % copy attributes from original file
    [~,~,~,natts] = netcdf.inqVar(ncid_in, varid_in);
    for a = 0:natts-1 % zero-offset
        attname = netcdf.inqAttName(ncid_in, varid_in, a);
        attrvalue = netcdf.getAtt(ncid_in, varid_in, attname);
        if strcmp(attname, '_FillValue')
            netcdf.defVarFill(ncid_out, varid_woa{v}, false, attrvalue);
        else
            netcdf.putAtt(ncid_out, varid_woa{v}, attname, attrvalue);
        end
    end
    netcdf.close(ncid_in)
end
netcdf.endDef(ncid_out)

% output of WOA09 depth levels
netcdf.putVar(ncid_out, varid_z, woa09_depths)

% determine indices of depth levels that will be kept
is_woa09_grid = ismember(woa18_depths, woa09_depths);

% copy values at selected depths for main variables
for v = 1:length(WOA_var_out)
    this_var = ncread(fn_woa_combined_var{v}, var_name_in{v});
    this_var_reduced = this_var(:,:,is_woa09_grid,:);
    netcdf.putVar(ncid_out, varid_woa{v},this_var_reduced)
end

netcdf.close(ncid_out)