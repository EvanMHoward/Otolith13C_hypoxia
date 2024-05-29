function [ WOAsub ] = subset_WOA( WOA, vnames, zmax )


%Author: Curtis Deutsch (cdeutsch@princeton.edu)
%Original version: 2021
%Current version: 2021

%Description: Select ony depths and variables of interest from gridded WOA
%fields.

%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 

v2struct(WOA)
clear WOA

% subset dimensions
Iz=find(z<=zmax);
z=z(Iz); dz=dz(Iz);
ztop=ztop(Iz);
zbot=zbot(Iz);
zedge=zedge(1:length(Iz)+1);

% subset tracer fields
for i=1:length(vnames)
    eval(['var=' vnames{i} ';']);
    var=var(:,:,Iz,:);
    eval([vnames{i} '=var;']);
end
if exist('po2','var')
    po2=po2(:,:,Iz,:);
end
clear var i 

% other fields
v3d=v3d(:,:,Iz);

WOAsub=v2struct;

end

