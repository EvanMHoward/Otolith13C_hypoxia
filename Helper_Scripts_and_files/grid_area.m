function [dx,dy,gridarea] = grid_area(lon,lat)
% GRID_AREA This function computes the grid area for each cell on a 
% spherical Earth.
% Author: unknown
% Original version: unknown
% Current version: August 27, 2024

r = 6371e3; % radius of earth in meters

dlat=diff(lat)*pi/180;
dlat(length(lat))=dlat(1);
dlon=diff(lon)*pi/180;
dlon(length(lon))=dlon(1);
dxfactor=cos(abs(lat*pi/180));

dx = nan(length(lat),length(lon));
dy = nan(length(lat),length(lon));

for j=1:length(lon)
    for i=1:length(lat)
        dx(i,j)=r*dlon(j)*dxfactor(i);
        dy(i,j)=r*dlat(i);
    end
end

dx=dx';
dy=dy';
gridarea=dx.*dy;
