%Author: Evan Howard (ehoward2@uw.edu)
%Original version: March 2024
%Current version: March 2024

%Description: This script fits an empirical model of d13CDIC following the
%functional form of Eide et al. 2017. 

%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 

%% Ocean fields
addpath('.../Helper_Scripts_and_files/'); %Location of helper scripts

zmax=200;
vnames={'temp','salt','o2','po4','aou'};
%World ocean atlas 2018
    load WOA
%Global ocean data assimilation project V2 (2016 climatology)
    load GLO
%WOA gridding
    load Goc

%Climatological AOU
o2sat=O2sol(WOA.salt,WOA.temp);
WOA.aou=o2sat-WOA.o2;

% subset Ocean for upper 200 m
WOAsub=subset_WOA(WOA,vnames,zmax);
Iz=WOAsub.Iz;

%% Generate d13DIC model
load('GLODAPv2.2023_Merged_Master_File.mat'); %Full GLodapV2 2023 release
I=G2depth<2e4; % maximum depth (include all depths)
X=[G2temperature, G2salinity, G2phosphate, G2aou];
[bDIC,seDIC,pDIC,modDIC,statsDIC] = stepwisefit(X(I,:),G2c13(I));
b=[bDIC;statsDIC.intercept];
dic13=b(1)*G2temperature(I)+b(2)*G2salinity(I)+b(3)*G2phosphate(I)+b(4)*G2aou(I)+b(5);
dif13=dic13-G2c13(I);
fmo=fitlm(dic13,G2c13(I)); %regress data to model (flipped x-y to force slope fit and evaluate uncertainty on 0 intercept)
seInt=fmo.Coefficients{1,2};seDIC(5,1)=seInt;sdDIC=seDIC;
%Matlab built-ins for fitlm and stepwiselm mislabeled: these uncertainties are the sqrt of covariance
%matrix, thus are standard deviations, not standard errors
mdl13DIC={b,sdDIC,pDIC,statsDIC};
save mdl13DIC mdl13DIC;

%% Plots to check output
figure;plot(G2c13(I),dic13,'k.');axis([-1 2 -1 2]); %predicted vs observed d13DIC

%Specific depths
wdic13=b(1)*WOA.temp+b(2)*WOA.salt+b(3)*WOA.po4+b(4)*WOA.aou+b(5);
pz=[0 100 200];
z=double(WOA.z);
x=G2longitude;y=G2latitude;
x(x<0)=x(x<0)+360;
cl=[-0.5 2];
figure(1); clf reset
for i=1:length(pz)
subplot(3,2,2*i-1)
I=abs(G2depth-pz(i))<20;
scatter(x(I),y(I),y(I)*0+50,G2c13(I),'filled')
xlim([0 360]);ylim([-80 80])
colorbar
clim(cl)
title(['Obs: z=' num2str(pz(i))])
subplot(3,2,2*i)
pcolor(WOA.x,WOA.y,nanmean(wdic13(:,:,find(z==pz(i)),:),4)');
shading flat;colorbar
title(['Mod: z=' num2str(pz(i))])
xlim([0 360]);ylim([-80 80])
clim(cl)
end



