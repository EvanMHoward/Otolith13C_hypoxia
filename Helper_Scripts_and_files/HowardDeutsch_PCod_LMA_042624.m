%Author: Evan Howard (ehoward2@uw.edu)
%Original version: January 2021
%Current version: October 2021

%Description: This script generates length-mass-age fits to the Pacific cod
%datasets after pairing with NOAA AFSC haul data.

%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 

%% Prepare workspace
clear all;close all;clc;
addpath('.../Helper_Scripts_and_files/'); %Location of helper scripts

swLMA=0; %if swLMA==1, fit function parameters for length-mass-age relationships
%% Load haul and fish size data
load HaulInfo; %Pacific cod sampled for otoliths
%Cruise,Haul,Specimen,Sex(unk),Length(mm),Weight(g),T_bottom,T_surface
l=HaulInfo(:,5);w=HaulInfo(:,6); %fish length (mm) and mass (g)

%Ages (ring counting, peak counting) [Kastelle et al. 2017]
age=nan(size(HaulInfo,1),2); %Read of manually, some differences between 18O and ring count age
age=[2,2;2,2;2,2;2,2;2,3;
    3,2;3,2;3,3;3,3;4,3;
    3,3;3,3;3,3;4,4;4,5;
    4,5;5,5;5,nan;2,2;3,3;
    4,4;3,3;4,3;3,3;5,5;
    4,5;4,4;4,nan;5,4;4,5;
    5,4;5,3;2,2;2,3;2,2;
    2,2;5,5;5,4;5,4;5,5];
age(isnan(age(:,2)),2)=age(isnan(age(:,2)),1); %duplicate or remove missing ages

%% Length-mass-age relationships
if swLMA==1
%fit Von Bertalanffy growth curves
funct_A=@(x,p)p(1).*(x.^p(2));pin_A=[1.7e-6 3.3];
[~,param_A,~,~,~,covparam_A,covresid_A,~,~,r2_A]=nlleasqr(l,w,pin_A,funct_A);
funct_B=@(x,p)p(1).*(1-exp(-p(2).*(x-p(3))));pin_B=[850 0.27 0.34];
[~,param_B1,~,~,~,covparam_B1,covresid_B1,~,~,r2_B1]=nlleasqr(age(:,1),l,pin_B,funct_B);
[~,param_B2,~,~,~,covparam_B2,covresid_B2,~,~,r2_B2]=nlleasqr(age(:,2),l,pin_B,funct_B);
funct_C1=@(x,p)p(1).*((1-exp(-param_B1(2).*(x-param_B1(3)))).^param_A(2));pin_C1=[7860];
funct_C2=@(x,p)p(1).*((1-exp(-param_B2(2).*(x-param_B2(3)))).^param_A(2));pin_C2=[3590];
[~,param_C1,~,~,~,covparam_C1,covresid_C1,~,~,r2_C1]=nlleasqr(age(:,1),w,pin_C1,funct_C1);
[~,param_C2,~,~,~,covparam_C2,covresid_C2,~,~,r2_C2]=nlleasqr(age(:,2),w,pin_C2,funct_C2);

lrng=[10:10:800];trng=[0.1:0.1:7];
subplot(131);plot(l,w,'.k');hold on;plot(lrng,funct_A(lrng,param_A),'-k');
ylabel('mass (g)');xlabel('length (mm)');
subplot(132);plot(age(:,1),l,'.k');hold on; plot(age(:,2),l,'.r');
plot(trng,funct_B(trng,param_B1),'-k');plot(trng,funct_B(trng,param_B2),'-r');
ylabel('length (mm)');xlabel('age (yrs)');
subplot(133);plot(age(:,1),w,'.k');hold on; plot(age(:,2),w,'.r');
plot(trng,funct_C1(trng,param_C1),'-k');plot(trng,funct_C2(trng,param_C2),'-r');
ylabel('mass (g)');xlabel('age (yrs; black is rings, red is peaks)');

LMAparams=cell(5,4);
%model A,B1,B2,C1,C2 by [params, parameter covariance (uncertainty), 
%residual covarance (error), and fit r^2].

LMAparams{1,1}=param_A;LMAparams{1,2}=covparam_A;LMAparams{1,3}=covresid_A;LMAparams{1,4}=r2_A;
LMAparams{2,1}=param_B1;LMAparams{2,2}=covparam_B1;LMAparams{2,3}=covresid_B1;LMAparams{2,4}=r2_B1;
LMAparams{3,1}=param_B2;LMAparams{3,2}=covparam_B2;LMAparams{3,3}=covresid_B2;LMAparams{3,4}=r2_B2;
LMAparams{4,1}=param_C1;LMAparams{4,2}=covparam_C1;LMAparams{4,3}=covresid_C1;LMAparams{4,4}=r2_C1;
LMAparams{5,1}=param_C2;LMAparams{5,2}=covparam_C2;LMAparams{5,3}=covresid_C2;LMAparams{5,4}=r2_C2;

save('LMAparams','LMAparams');
save('LMAfunctions','funct_A','funct_B','funct_C1','funct_C2');

clearvars lrng trng param_A param_B1 param_B2 param_C1 param_C2;
clearvars covparam_A covparam_B1 covparam_B2 covparam_C1 covparam_C2;
clearvars covresid_A covresid_B1 covresid_B2 covresid_C1 covresid_C2;
clearvars r2_A r2_B1 r2_B2 r2_C1 r2_C2 pin_A pin_B pin_C1 pin_C2;

else
    load LMAparams;
    load LMAfunctions;
end
%Looks to me like the 18O-peak based age model has some length-mass
%outliers that are unlikely. Noise in the measurements may be obscuring 
% true peaks.
