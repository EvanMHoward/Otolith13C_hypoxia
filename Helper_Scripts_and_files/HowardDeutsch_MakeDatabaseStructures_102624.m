%Author: Evan Howard (ehoward2@uw.edu)
%Original version: February 2022
%Current version: October 2024

%Description: The other scripts call structures 'WOA', 'GLO', and
%'TrophicLevel'. The glodap field is regenerated from the original 
% (online) and the other two are regenerated from more compact structures.

%Sources of data:
%TrophicLevel: %These data were extracted from the FishBase 'Ecology' table 
%in RFishBase (downloaded May 2022). No support is provided here for R
%functions, instead the relevant variables are stored with the provided
%scripts in this repository.

%GLO: GLODAPv2 2016 climatologies ('GLODAPv2.2016b_MappedClimatologies.tar.gz
%downloaded from https://glodap.info/index.php/mapped-data-product/).

%WOA: WOA2018 climatologies hosted by NCEI, 
% https://www.ncei.noaa.gov/data/oceans/woa (see
% 'Driver_HowardEtAl_otoliths_WOA.m')

%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 

%% Prepare workspace
clear all;close all;clc;
addpath(genpath('.../Helper_Scripts_and_files')); %Location of helper scripts
                                                %modify as needed

gdir='~/GLODAPv2.2016b_MappedClimatologies/'; %Modify to the location user 
%                       has stored unzipped GLODAP climatology netcdf files

%% Regenerate the FishBase derived 'TrophicLevel' structure
%The default value for this application is the Foodweb based trophic level 
%estimate ('FoodTroph'), with missing values filled in using the 
%dietary-study (literature-based) estimates ('DietTroph') when FoodTroph is 
%missing but DietTroph is present. The former is a more complete list, but 
%overall the two distributions are indistinguishable.

if ~isfile("TrophicLevel.mat")
load TrophicLevelTable;
a=TrophicLevelTable{:,1};
b=TrophicLevelTable{:,2};
c=TrophicLevelTable{:,3};
TrophicLevel={a,b,c};
save TrophicLevel.mat TrophicLevel -v7.3
end

%% Regenerate Glodap 2016 climatology reshaped to WOA gridding ('GLO')

if ~isfile("GLO.mat")
    ddir=gdir;
    fname='GLODAPv2.2016b.temperature.nc';
    Temp=ncread([ddir fname],'temperature');
    fname='GLODAPv2.2016b.salinity.nc';
    Salt=ncread([ddir fname],'salinity');
    fname='GLODAPv2.2016b.TCO2.nc';
    DIC=ncread([ddir fname],'TCO2');
    fname='GLODAPv2.2016b.TAlk.nc';
    Alk=ncread([ddir fname],'TAlk');
    fname='GLODAPv2.2016b.silicate.nc';
    Si=ncread([ddir fname],'silicate');
    fname='GLODAPv2.2016b.PO4.nc';
    PO4=ncread([ddir fname],'PO4');
    
    % pressure
    y=ncread([ddir fname],'lat');
    x=ncread([ddir fname],'lon');
    z=ncread([ddir fname],'Depth');
    [gx,gy,gz]=ndgrid(x,y,z);
    Prs=sw_pres(gz,gy);
    clear g*;
    
    xsh=20;
    x=circshift(x,xsh);
    Temp=circshift(Temp,xsh,1);
    Salt=circshift(Salt,xsh,1);
    DIC=circshift(DIC,xsh,1);
    Alk=circshift(Alk,xsh,1);
    
    A=CO2SYS(Alk(:),DIC(:),1,2,Salt(:),Temp(:),nan,Prs(:),nan,Si(:),PO4(:),0,0,1,15,3,2,2);
    A(A==-999)=nan;
    
    pco2 = reshape(A(:,4),size(Temp));
    pH = reshape(A(:,3),size(Temp)); 
    OmA = reshape(A(:,18),size(Temp)); 
    OmC = reshape(A(:,17),size(Temp)); 
    
    clear A ddir fname xsh i j k;
    
    GLO=v2struct;
    save GLO.mat GLO -v7.3
end

%% Regenerate WOA 2018 climatology structure ('WOA')
%This section calls a driver and associated helper scripts, which directly
%download and process the relevant WOA data. See
%'Driver_HowardEtAl_otoliths_WOA.m' for additional details and to change
%default filepaths.
%NOTE: WOA concentrations are reported in mmol m-3, but this script
%converts all chemical units to umol kg-1 for compatability with other
%scripts used in this analysis.
clear all; clc;
run('Driver_HowardEtAl_otoliths_WOA.m');
