%Author: Evan Howard (ehoward2@uw.edu)
%Original version: February 2022
%Current version: September 2022

%Description: This scripts pairs the trophic level and d13Com-d13Cdic(sw)
%for each of the species considered in the global analysis.

%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 

%% Prepare workspace
clear all;close all;clc;
addpath('.../Helper_Scripts_and_files/'); %Location of helper scripts
load InterSpecies_b;load TrophicLevel; %Trophic level database from FishBase, downloaded September 2022

%Find offset by species
for i=1:length(C_IS)
aa=C_IS{i,1}(:,5)-C_IS{i,1}(:,4);bb=nanmean(aa);cc=nanstd(aa);
if cc<1e-10;cc=nan;end
D_OMminDIC(i,:)=[bb,cc];
end

%pair to trophic level
FoodTroph=TrophicLevel{:,2};
i=isnan(FoodTroph);
FoodTroph(i)=TrophicLevel{1,3}(i);
spec2=string(TrophicLevel{:,1});

clearvars TL aa;
for i=1:length(C_IS)
    spec1=string(C_IS{i,2});
    j=strcmp(spec1,spec2);
    if sum(j)>0
        aa=nanmean(FoodTroph(j));
        TL(i,1)=aa;
    else TL(i,1)=NaN;
    end
    if i==14;TL(i,1)=4.5;end %looked up manually on online FishBase, missing in their export table
end

clearvars aa;
aa=[TL,D_OMminDIC];
%only 5,6,8 are indpendent cod datsets, keeping these separate for now but
%removing duplicate value from 7.
aa=aa([1:6,8:end],:); 
%includes two estimated trophic levels without associated enrichments:
%could exclude those by ending at index 12, but only changes trophic level
%from 4.08 to 4.10
x=aa(:,1);y=aa(:,2);w=1./(aa(:,3).^2);
[nanmean(x),nanstd(x),nanmean(y),nanstd(y)]
