%Author: Evan Howard (ehoward2@uw.edu)
%Original version: December 2020
%Current version: April 2024

%Description: This script compares temperature-dependent hypoxia based 
% models to observed Atlantic cod otolith data, as well as an analysis of 
% how the uncertainty in the metabolic to seawater carbon ratio changes 
% with respect to measured isotopic compositions.

%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 

%% Prepare workspace
clear all; close all; clc;
addpath('.../Helper_Scripts_and_files/'); %Location of helper scripts
load InterSpecies_b

%Default parameters
%endmembers
e_om=1; %approximate fractionation of blood DIC versus otolith organic matter
e_om2=e_om; %blood DIC vs diet/tissue--may be different than otolith OM
%Trophic enrichment of 1.5 per mil (1-2 per mil typical for marine fish) 
% reported by Sweeting et al. 2007.

%Set del_dic and del_om initially, but this can be changed as needed with 
% different datasets. 
%Here considering Chung et al. 2019 reared cod dataset: other reared and
%wild cod datasets generate similar results, but require normalization for
%different metabolic endmembers in order share a single model.
del_dic=mean(C_IS{5}(:,4));del_om=mean(C_IS{5}(:,5))+e_om;

%environmental chemistry
pCO2_meas=370e-6; %approximate annual mean upper ocean pCO2
pO2_meas=0.209; %approximate annual mean upper ocean pO2
spbC='CO2aq'; %blood carbon species for mass balance

%Van't Hoff equation constants
kb=8.617e-5;T0=273.15;Tref=15;

%Atlantic cod hypoxia parameters from prior work (Deutsch et al. 2015)
Pc=2.4; %Atlantic cod estimate 
del=0.82; %demand (metabolism) only allometric exponent
eps=-0.06; %supply-demand allometric exponent (note direction convention)
%%%NOTE: if exponents delta and sigma are not biomass normalized, then the
%%%ratio of demand to supply for biomass normalized rates is
%%%B^(delta-1)/B^(sigma-1)=B^(delta-sigma)=B^-eps. i.e., the same result as
%%%as when rates are not biomass normalized.
Em=0.519; %demand (metabolism) only temperature exponent
Eo=0.347; %demand-supply temperature exponent (note different direction convention)
%convert this to Eo at Tref
dEdT=0.01; %temperature dependence of Eo (Atlantic cod)
ET=5; %temperature at which standard Eo is observed
Eo=Eo+dEdT*(Tref-ET); % Eo value at Tref
Eor=Eo; %saved reference
Vh=0.0621; %hypoxia vulnerability (atm) at 15 degrees (C)
Bref=250; %Median mass of all cod in from these four datasets (2 reared and 2 wild)
%Allometric reference (Plante et al. 1998) has median mass around 1200g.,
%however, here we are trying to look at variation with respect to the
%datasets from which otoliths were taken, so lower mass may be appropriate.
%Median of plotted data is about 250g.

%ranges of temperature and mass to consider          
trange=0:0.5:30; %Temperature range considered, leaving salinity fixed
brange=logspace(0,4,100); %Mass range considered, 1g to 10 kg 
%(have to include larval sizes for wild cod data)
brangeN=brange./Bref; %Normalized mass compared to reference

%switches
swF=1; %use different blood chemistry than seawater chemistry
swStoi=1; %add factors to qd or qs 
        fd=0.5; %Vh_old*0.5, adjusts Vh from assumed pO2_blood=0 at end of 
        %circulation to average pO2_blood=0.5*pO2_environment (to compare
        %box model with prior respirometry assumptions)
        fs=1; %factor of CO2 ventilation related to metabolism
swEoT=1; %Switch to include approximated temperature sensitivity of Eo based on 
        %transition from ventilation limited to diffusion limited supply
%% Aqueous chemistry and switches
%Arrhenius temperature term
TT=(1/kb).*(1./(trange+T0)-1./(Tref+T0));

%Environmental gas partial pressures
pCO2ref=pCO2_meas;pO2ref=pO2_meas;

%Assumed pH and salinity/ionic strength of environment and fish fluids
%This is a minor effect so don't need to be very accurate on pH values
pH_oc=8.07;TpHref_oc=10; %Reference pH (total scale) and temperature, based on GLODAP
pH_blood=7.85-0.15;TpHref_blood=TpHref_oc; %pH_NBSscale of fish blood as 
%low as 7.3, but generally around 7.7-8.0; -0.15 approximates scale difference
S_oc=34;S_f=9; %approximate ocean salinity, as well as S with roughly 
% equivalent ionic stregth of fish plasma or endolymph
if swF~=1
    S_f=S_oc;pH_blood=pH_oc; %overwrite fish blood with seawater chemistry
end

%Temperature sensitivity of pH
pHe=pH_T(trange,S_oc,pH_oc,TpHref_oc);pHe=pHe';
pHf=pH_T(trange,S_f,pH_blood,TpHref_blood);pHf=pHf';

%Carbon dioxide solubility
[K0e K1e K2e]=carbeq(trange,S_oc);K0e=K0e';K1e=K1e';K2e=K2e'; %In seawater
[K0f K1f K2f]=carbeq(trange,S_f);K0f=K0f';K1f=K1f';K2f=K2f'; %In fish bodily fluids

%Carbonate system partioning
He=10.^(-1.*pHe);De=(He.^2)+He.*K1e+K1e.*K2e;
X0e=(He.^2)./De;X1e=He.*K1e./De;%X2=1-X0-X1;
X2e=K1e.*K2e./De;
Hf=10.^(-1.*pHf);Df=(Hf.^2)+Hf.*K1f+K1f.*K2f;
X0f=(Hf.^2)./Df;X1f=Hf.*K1f./Df;%X2=1-X0-X1;
X2f=K1f.*K2f./Df;

%Oxygen solubility
Khe=O2sol(S_oc,trange)./(pO2ref).*1e-6;Khe=Khe'; %mol kg-1 atm-1 of O2
Khf=O2sol(S_f,trange)./(pO2ref).*1e-6;Khf=Khf'; 

%Diffusivities in water (all in 10-4 m2 s-1)
[Dco2e,Dhco3e,Dco3e]=co2_diffusion(S_oc,trange);Dco2e=Dco2e';Dhco3e=Dhco3e';Dco3e=Dco3e';
[Dco2f,Dhco3f,Dco3f]=co2_diffusion(S_f,trange);Dco2f=Dco2f';Dhco3f=Dhco3f';Dco3f=Dco3f';
[~,~,~,~,Do2e,~,~]=gas_diffusion(S_oc,trange);Do2e=Do2e';
[~,~,~,~,Do2f,~,~]=gas_diffusion(S_f,trange); Do2f=Do2f';

%Isotopic endmembers (for reared cod data)
stn=0.011180;
def=del_dic;dmf=del_om;
Ref=((def)./1000+1).*stn;Rmf=((dmf)./1000+1).*stn; %fixed carbon endmembers
tmpRe=nan(length(trange),1);
for i=1:length(trange)
[~,tReG,~,~,~,~]=C13speciation(Ref,pHe(i),K1e(i),K2e(i),trange(i),'DIC'); %CO2aq,e   
tmpRe(i)=tReG;clearvars tReG;
end
[~,iT]=min(abs(trange-Tref));
tmpRe=tmpRe';Reref=tmpRe(iT);

%Stoichiometry
qd=1;
qs=K0e.*Dco2e./(Khe.*Do2e);
if swStoi==1
qs=qs.*fs;qd=qd.*fd;
end

%Temperature dependence of Eo
if swEoT==1
Eo=Eor+dEdT.*(trange-Tref);
end
%% Temperature sensitivity
%Plot as change versus reference instead of absolute, to make presentation
%more consistent across panels. 
[~,iT]=min(abs(trange-Tref)); %intercept of Ed and Eo to use
Reref=tmpRe(iT);
[~,iT2]=min(abs(trange-8.5)); %centered data

% Reared Atlantic cod
clearvars RoT;
RoT=nan(2,size(trange,2));
qsT=qs';
Pctmp=1; %reared cod appear to have PhiCrit close to 1 (as evaluated later)
%Assuming a constant here lets us focus on the temeprature driven trend
%only, cant solve for every parameter at once.

%full model
RoCO2=(Rmf.*Vh.*Pctmp./pCO2ref.*qd./qsT.*exp(-Eo.*TT)+tmpRe)./(Vh.*Pctmp./pCO2ref.*qd./qsT.*exp(-Eo.*TT)+1);
for i=1:length(trange)
[~,~,~,~,R13Arag,~]=C13speciation(RoCO2(i),pHf(i),K1f(i),K2f(i),trange(i),spbC);
RoT(1,i)=R13Arag;
end
clearvars R13Arag RoCO2;
%Metabolism dependence only
RoCO2=(Rmf.*Vh.*Pctmp./pCO2ref.*qd./qsT.*exp(-Em.*TT)+tmpRe)./(Vh.*Pctmp./pCO2ref.*qd./qsT.*exp(-Em.*TT)+1);
for i=1:length(trange)
[~,~,~,~,R13Arag,~]=C13speciation(RoCO2(i),pHf(i),K1f(i),K2f(i),trange(i),spbC);
RoT(2,i)=R13Arag;
end
clearvars qsT R13Arag RoCO2;

%Isotopic fractionations only with constant carbonate system (but including
%the change in Rw from temperature dependent fractionation relative to DIC)
qsT=qs(iT);
RoCO2=(Rmf.*Vh.*Pctmp./pCO2ref.*qd./qsT+tmpRe)./(Vh.*Pctmp./pCO2ref.*qd./qsT+1);
[~,~,~,~,R13Arag,~]=C13speciation(RoCO2,pHf(iT),K1f(iT),K2f(iT),trange,spbC);
RoT(3,:)=R13Arag;clearvars qsT R13Arag RoCO2;
%Stoichiometry (Qs) changes with temperature only
qsT=qs';
RoCO2=(Rmf.*Vh.*Pctmp./pCO2ref.*qd./qsT+Reref)./(Vh.*Pctmp./pCO2ref.*qd./qsT+1);
[~,~,~,~,R13Arag,~]=C13speciation(RoCO2,pHf(iT),K1f(iT),K2f(iT),trange(iT),spbC);
RoT(4,:)=R13Arag;clearvars qsT R13Arag RoCO2;

%find offsets
aa=(RoT(1,:)./stn-1).*1000;
%save offset for normalization and recenter plot so 0 change is near 8.5C
%for this dataset (but intercept between Ed and Eo still at 15C to
%illustrate slope divergence).
offst1=aa(iT); %15C centered along Eo line
offst2=aa(iT2); %8.5C centered along Eo line

%Plot
subplot(221);
%plot chemistry only terms
%Isotopic fractionations including change in Rw
aa=(RoT(3,:)./stn-1).*1000;aa3=aa-aa(iT)-(offst2-offst1); 
patch([trange,fliplr(trange)],[ones(1,length(aa3)).*(aa3(iT)),fliplr(aa3)],'blue','EdgeColor','b');hold on;
%Qs (pH turns out to be tiny in net)
aa=(RoT(4,:)./stn-1).*1000;cl=aa-aa(iT);aa4=aa-aa(iT)-(offst2-offst1);
fl=aa3;
patch([trange,fliplr(trange)],[fl,fliplr(cl+fl)],'red','EdgeColor','r');hold on;

%plot full model with either Eo or Ed as metabolic temperature sensitivity
aa=(RoT(1,:)./stn-1).*1000;
aa=aa-aa(iT)-(offst2-offst1); %normalize to 15C for Ed & Eo intercept, 
%then offset entire line extra factor so everything moves in parallel to
%center data on 0.
aa1=aa;
aa=(RoT(2,:)./stn-1).*1000;
aa=aa-aa(iT)-(offst2-offst1);  %normalize to 15C for Ed & Eo intercept, 
%then offset entire line extra factor so everything moves in parallel to
%center data on 0.
aa2=aa;
plot(trange,aa2,'--k','LineWidth',2,'Color',[0.5 0.5 0.4]);hold on; %Em
plot(trange,aa1,'-k','LineWidth',2);hold on; %Eo
clearvars i i1;
xx=C_IS{5}(:,2);yy=C_IS{5}(:,1);
i1=findgroups(xx);
for i=1:length(unique(i1))
    x=mean(xx(i1==i));
    y=mean(yy(i1==i));
    ypos=max(yy(i1==i))-y;
    yneg=y-min(yy(i1==i));
    y=y-(offst2);
    errorbar(x,y,yneg,ypos,'.k','MarkerSize',16,...
        'MarkerEdgeColor','k','CapSize',0,'LineWidth',1.5);
    clearvars x y ypos yneg;
end
clearvars i i1;
xlabel('Temperature (^{o}C)');ylabel('\Delta\delta^{13}C_{o}');
xlim([0 20]);ylim([-2 2]);
set(gca,'Layer','top');

% Check stats on the two slopes
%aa1 is model with Eo, aa2 is model with Ed, xx is data T, yy-ofsfst2 is
%data 13C

%Model values at observed temperatures
clearvars it iit;
it=nan(length(xx),1);
for i=1:length(xx)
    [~,iit]=min(abs(xx(i)-trange));
    it(i)=iit;
end

yd1=((yy-offst2)-aa1(it)');yd2=((yy-offst2)-aa2(it)'); %data difference around models is close enough to normal for this purpose
s1=std(yd1);s2=std(yd2); %standard deviation of each
L1=(-length(xx)./2.*log(2.*pi())-length(xx).*log(s1)-1./(2.*(s1.^2)).*sum(yd1.^2)); %log-likelihood
L2=(-length(xx)./2.*log(2.*pi())-length(xx).*log(s2)-1./(2.*(s2.^2)).*sum(yd2.^2)); 
q1=2;q2=1; %#number of parameters (hypoxia traits, other inputs treated as constants); 
% %only difference matters, and main model has dEo/dT
AIC1=-2*L1+2*q1+2*q1*(q1+1)/(length(xx)-q1-1);
AIC2=-2*L2+2*q2+2*q2*(q2+1)/(length(xx)-q2-1); %%Akaike information criteria, with finite sample correction
am=min(AIC1,AIC2);B1=exp((am-AIC1)/2);B2=exp((am-AIC2)/2); %Relative difference from minimum
pAIC1=B1/(B1+B2); pAIC2=B2/(B1+B2); %Akaike weights (conditional probability)
%Even with the extra dEo/dT term, Eo model is vastly superior to Ed
%What is equivalent p value?
lambda=(AIC2-AIC1);
p = chi2cdf(lambda,(q1-q2),'upper'); %p value for X2 distribution
%% Biomass dependence
%Generate model prediction for precise endmembers of each sample
clearvars tRoT;tRoT=cell(4,2);
for j=1:4 %reared and wild datasets
    if j==1 id=5;Pctmp=1; %Chung reared
    elseif j==2 id=7;Pctmp=1; %Gao reared
    elseif j==3 id=6;Pctmp=2.4; %Jamieson wild %reset Pctmp for wild cod later when Panel 1 finished
    elseif j==4 id=8;Pctmp=2.4; %Gao wild
    end
    d_dic=C_IS{id}(:,4);dmeas=C_IS{id}(:,1);
    t=C_IS{id}(:,2);b=C_IS{id}(:,3);
    if j==1 | j==2
    d_om=C_IS{id}(:,5)+e_om;
    else
    d_om=C_IS{id}(:,5)+e_om2;
    end
    
    %Do temperature model for each endmember, but same standard trange
    for k=1:length(t)
    tRef=((d_dic(k))./1000+1).*stn;tRmf=((d_om(k))./1000+1).*stn; %datum specific carbon endmembers
    ttmpRe=nan(length(trange),1);
    for i=1:length(trange)
    [~,ReG,~,~,~,~]=C13speciation(tRef,pHe(i),K1e(i),K2e(i),trange(i),'DIC'); %CO2aq,e   
    ttmpRe(i)=ReG;clearvars ReG;
    end
    ttmpRe=ttmpRe';  
    qsT=qs';
    RoCO2=(tRmf.*Vh.*Pctmp./pCO2ref.*qd./qsT.*exp(-Eo.*TT)+ttmpRe)./...
        (Vh.*Pctmp./pCO2ref.*qd./qsT.*exp(-Eo.*TT)+1);
    for i=1:length(trange)
    [~,~,~,~,R13Arag,~]=C13speciation(RoCO2(i),pHf(i),K1f(i),K2f(i),trange(i),spbC);
    tRoT{j,1}(k,i)=R13Arag;
    end 
    
    %What is the anomaly of the measured data versus the temperature-dependent
    %model for each datum's specific endmembers?
    iTk=find(trange==t(k)); %index of measured temperature
    dmod=(tRoT{j}(k,iTk)./stn-1).*1000; %del13Coto modeled for temperature dependence only
    tRoT{j,2}(k,1)=dmeas(k)-dmod; %this is anomaly between measured and modeled
    tRoT{j,2}(k,2)=dmeas(k);
    tRoT{j,2}(k,3)=t(k);
    tRoT{j,2}(k,4)=b(k);
    
    clearvars R13Arag RoCO2 dmod iTk tRef tRmf ttmpRe;
    clearvars dmod dcorr;
    end %length of data vector in each dataset
    clearvars  t b d_dic d_dom d_meas;
end %dataset

%Check a linear fit of the available data
clearvars o b;
o=cat(1,tRoT{1,2}(:,1),tRoT{2,2}(:,1),tRoT{3,2}(:,1),tRoT{4,2}(:,1));
b=log10(cat(1,tRoT{1,2}(:,4),tRoT{2,2}(:,4),tRoT{3,2}(:,4),tRoT{4,2}(:,4)));
mdlM=fit(b,o,'poly1','Robust','LAR');cima=predint(mdlM,log10(brange),0.95,'functional','on');
omdl=o;bmdl=b;
dfm=2-1;dfe=length(omdl)-2;SSM=sum((mdlM(bmdl)-mean(omdl)).^2);SSE=sum((omdl-mdlM(bmdl)).^2);
Fmdl=(SSM/dfm)/(SSE/dfe);pmdl=1-fcdf(Fmdl,dfm,dfe);
%using least absolute residuals weighting because otherwise varying sampling 
% density biases fit towards better sampled medium to high mass.

%Mass dependent model for mean of reared cod data
%Note, to use with anomaly plot, need to compare to model that has both
%temp and mass
RoB=nan(2,length(brangeN));
qsT=qs(iT);
Pctmp=1; %reared cod appear to have PhiCrit close to 1 (evaluated later)
%Assuming a constant here lets us focus on the mass driven trend
%only, cant solve for every parameter at once.

%supply-demand
RoCO2=(Rmf.*Vh.*Pctmp./pCO2ref.*qd./qsT.*(brangeN.^-eps)+Reref)./...
    (Vh.*Pctmp./pCO2ref.*qd./qsT.*(brangeN.^-eps)+1);
[~,~,~,~,R13Arag,~]=C13speciation(RoCO2,pHf(iT),K1f(iT),K2f(iT),trange(iT),spbC);
RoB(1,:)=R13Arag;clearvars R13Arag RoCO2;
clearvars R13Arag RoCO2;
%demand only
RoCO2=(Rmf.*Vh.*Pctmp./pCO2ref.*qd./qsT.*(brangeN.^(del-1))+Reref)./...
    (Vh.*Pctmp./pCO2ref.*qd./qsT.*(brangeN.^(del-1))+1);
[~,~,~,~,R13Arag,~]=C13speciation(RoCO2,pHf(iT),K1f(iT),K2f(iT),trange(iT),spbC);
RoB(2,:)=R13Arag;clearvars R13Arag RoCO2;
clearvars R13Arag RoCO2;

load Allometry.mat; %Bg: centered log10 mass in g, delta, sigma, 
% with epsilon=delta-sigma (opposite convention used in this work,
% sigma-delta)
B10c=abs(diff(Bg)./2)+Bg(2:end);clearvars Bg sigma delta;
B10c=B10c(72:end)';eps_sc=-1.*epsilon(72:end)';clearvars epsilon;
BcN=(10.^B10c)./Bref;
RoCO2=(Rmf.*Vh.*Pctmp./pCO2ref.*qd./qsT.*(BcN.^-eps_sc)+Reref)./...
    (Vh.*Pctmp./pCO2ref.*qd./qsT.*(BcN.^-eps_sc)+1);
[~,~,~,~,R13Arag,~]=C13speciation(RoCO2,pHf(iT),K1f(iT),K2f(iT),trange(iT),spbC);
RoBcont=R13Arag;clearvars R13Arag RoCO2;

%Plot the reared cod model (mean endmembers), and the individual 
%data-model comparisons (specific endmembers for each datum)
subplot(222);
for jj=1:4
    if jj==1 clr='k';
    elseif jj==2 clr='k';
    elseif jj==3 clr='r';
    elseif jj==4 clr='r';
    end
    xx=log10(tRoT{jj,2}(:,4)); %log10 mass of fish
    yy=tRoT{jj,2}(:,1); %Anomaly of measured - temperature only model
    if jj==1 | jj==2 %reared
    plot(xx,yy,'o','MarkerSize',4,'MarkerEdgeColor',clr,'MarkerFaceColor','k');hold on;
    else
    plot(xx,yy,'o','MarkerSize',4,'MarkerEdgeColor',clr,'MarkerFaceColor','w');hold on;
    end
end

aa=(RoB(1,:)./stn-1).*1000-(RoT(1,iT)./stn-1).*1000; 
ir=(brange>=1) & (brange<=10000);
plot(log10(brange(ir)),aa(ir),'-k','LineWidth',2);hold on; %Supply and demand
aa=((RoB(2,:))./stn-1).*1000-(RoT(1,iT)./stn-1).*1000;
plot(log10(brange),aa,'--','LineWidth',2,'Color',[0.5 0.5 0.4]); %Demand only
xlabel('Mass (g)');ylabel('\Delta\delta^{13}C_{o}');
xlim([-0.1 4]);ylim([-2 2]);
set(gca,'Layer','top');

%Compare fit vs assumed model other than Eps=0
beps=log10(cat(1,tRoT{1,2}(:,4),tRoT{2,2}(:,4),tRoT{3,2}(:,4),tRoT{4,2}(:,4))); %now normalized mass
oeps=o-(-eps).*beps;
odel=o-(del-1).*beps;
mdlMeps=fit(beps,oeps,'poly1','Robust','LAR');
dfm=2-1;dfe=length(oeps)-2;SSM=sum((mdlMeps(bmdl)-mean(oeps)).^2);SSE=sum((oeps-mdlMeps(bmdl)).^2);
Fmdleps=(SSM/dfm)/(SSE/dfe);pmdleps=1-fcdf(Fmdleps,dfm,dfe);
mdlMdel=fit(beps,odel,'poly1','Robust','LAR');
dfm=2-1;dfe=length(odel)-2;SSM=sum((mdlMdel(bmdl)-mean(odel)).^2);SSE=sum((odel-mdlMdel(bmdl)).^2);
Fmdldel=(SSM/dfm)/(SSE/dfe);pmdldel=1-fcdf(Fmdldel,dfm,dfe);

%% Activity level
%Infer active to resting supply to demand ratios, PhiCrit, from otolith
%isotopic ratio. Leave off allometric dependence since available evidence
%is that if present it cannot be robustly distinguished from zero-dependence.
        
%similar to section above, but just calculate for mean of dataset types  
clearvars PcT;PcT=cell(2,2);clearvars kTi gTi;
Pcvec=[0.1:0.1:15]; %what are standard values of PhiCrit to evaluate

for j=1:2 %reared and wild datasets
    if j==1 id=[5,7];k_t=[3,1];%Chung and Gao reared
    elseif j==2 id=[6,8];k_t=[3,1];%Jamieson and Gao wild
    end
    clearvars d_dic d_om d_meas t b;
    for g=1:2 %temperature groups
        t=C_IS{id(g)}(:,2);
        [gT,kgT]=kmeans(t,k_t(g),'Distance','cityblock','Start','plus'); %temperature
        [sortd,sorto]=sort(kgT,'ascend');
        if g==1
            clearvars kTi gTi;
            for k=1:length(sorto)
                kTi(k,1)=sortd(k);
                gTi(gT==sorto(k))=k;
            end
            if size(gTi,1)==1
                gTi=gTi';
            end
        else
            clearvars tkTi tgTi;
            tkTi=t;tgTi=4.*ones(length(t),1); 
            kTi=cat(1,kTi,tkTi);
            gTi=cat(1,gTi,tgTi);
            clearvars tkTi tgTi;
        end
        clearvars sorto sortd  kgT gT k;
    end
    %append all datasets of same type
    d_dic=cat(1,C_IS{id(1)}(:,4),C_IS{id(2)}(:,4));
    if j==1
    d_om=cat(1,C_IS{id(1)}(:,5)+e_om,C_IS{id(2)}(:,5)+e_om);
    else
    d_om=cat(1,C_IS{id(1)}(:,5)+e_om2,C_IS{id(2)}(:,5)+e_om2);    
    end
    d_meas=cat(1,C_IS{id(1)}(:,1),C_IS{id(2)}(:,1));
    t=cat(1,C_IS{id(1)}(:,2),C_IS{id(2)}(:,2));
    b=cat(1,C_IS{id(1)}(:,3),C_IS{id(2)}(:,3));

    %Do temperature model for each endmember, but same standard trange
    %Similar to above section, but this time use measured Ro to solve for
    %PhiCrit: X=Vh*f(T)*f(B)/pCO2, and PhiCrit=(Re-Ro)/(Ro-Rm)*1/X where
    %Re and Ro are in terms of CO2aq and Rm in terms of DIC. Note, 
    %Cr/Ce=(Re-Ro)/(Ro-Rm) as well.
    for k=1:length(t)
    tRef=((d_dic(k))./1000+1).*stn;tRmf=((d_om(k))./1000+1).*stn; %datum specific carbon endmembers
    RoArag=((d_meas(k))./1000+1).*stn; %measured aragonite ratio    
    ttmpRe=nan(length(trange),1);
    ttmpRo=nan(length(trange),1);
    for i=1:length(trange)
    [~,ReG,~,~,~,~]=C13speciation(tRef,pHe(i),K1e(i),K2e(i),trange(i),'DIC'); %CO2aq,e   
    ttmpRe(i)=ReG;clearvars ReG;
    %Find Ro measured as CO2aq for any given temperature
    [~,RoCO2aq,~,~,~,~]=C13speciation(RoArag,pHf(i),K1f(i),K2f(i),trange(i),'Arag');
    ttmpRo(i)=RoCO2aq;clearvars RoCO2aq;
    end
    ttmpRe=ttmpRe'; %ttmpRe is Re (CO2aq) endmember, tRmf is Rm (DIC) endmember
    ttmpRo=ttmpRo'; %ttmpRo is Ro (CO2aq) measured       
    qsT=qs';
    %Calculate PhiCrit from other values
    iTk=find(trange==t(k)); %index of measured temperature (works so long as exact T is in trange)
    X=(Vh/pCO2ref*qd/qsT(iTk)*exp(-Eo(iTk)*TT(iTk)));
    %R can have numerical errors given very small values close
    %together in denominator (within measurement errors of either value)
    R=((ttmpRe(iTk)-ttmpRo(iTk))/(ttmpRo(iTk)-tRmf)); 
    Rdenom=(ttmpRo(iTk)./stn-1).*1000-(tRmf./stn-1).*1000; % As del difference
    Pc_m=R/X;
    PcT{j,1}(k,1)=Pc_m;
    PcT{j,1}(k,2)=d_meas(k);
    PcT{j,1}(k,3)=gTi(k);
    PcT{j,1}(k,4)=t(k);
    PcT{j,1}(k,5)=R; %Cr/Ce
    PcT{j,1}(k,6)=Rdenom; %delint-delmet
    PcT{j,1}(k,7)=d_om(k); %delmet (=delom+e_om)
    PcT{j,1}(k,8)=(ttmpRe(iTk)./stn-1).*1000; %delw
    PcT{j,1}(k,9)=(ttmpRo(iTk)./stn-1).*1000; %delint
    PcT{j,1}(k,13)=qd./qsT(iTk); %Q(T)
    PcT{j,1}(k,14)=Eo(iTk); %expected Eo with dEo/dT
    PcT{j,1}(k,15)=b(k); %biomass
    PcT{j,1}(k,16)=d_dic(k); %del13w(DIC)

    %What would relevant values be at normalized temperature and same d13Cw?
    Tnew=7;iTn=find(trange==Tnew);
    %What is ratio of PCmet/PCw at new temperature?
    % R(7)=R(T)*Phic(7)/Phic(T)*exp(...7)/exp(...T)*Q(7)/Q(T)
    %Could assume either that R/Phic is preserved, or that Phic is constant
    %(limited evidence for variation in PhiCrit with T)
    Rnew=R.*exp(-Eo(iTn)./kb.*(1./(Tnew+T0)-1./(Tref+T0)))./exp(-Eo(iTk)./kb.*(1./(t(k)+T0)-1./(Tref+T0)));
    %Use new assumed constant water composition
    Wnew=(0./1000+1).*stn; %mean water across datasets is ~+0.5 per mil (so assuming zero is reasonable)
    [~,ReG2,~,~,~,~]=C13speciation(Wnew,pHe(iTn),K1e(iTn),K2e(iTn),trange(iTn),'DIC'); %CO2aq,e  
    %And organic matter composition still varies as d_om(k) or tRmf
    %Then dCint= (R*dCmet+dCw)/(R+1), and can be converted to d13Coto
    RCint=(Rnew*tRmf+ReG2)/(Rnew+1); 
    [~,~,~,~,Roto2,~]=C13speciation(RCint,pHf(iTn),K1f(iTn),K2f(iTn),trange(iTn),'CO2aq');
    Doto2=(Roto2./stn-1).*1000;

    PcT{j,1}(k,10)=Doto2; %new otolith composition at constant T, varied d13met
    PcT{j,1}(k,11)=Tnew; %new temperature
    PcT{j,1}(k,12)=(ReG2./stn-1).*1000; %new shared water composition

    %Now do same excercise with measured d13w, but constant d13met (-19 per mil)
    dm2=-19;Rm2=(dm2/1000+1)*stn;
    [~,ReG3,~,~,~,~]=C13speciation(tRef,pHe(iTn),K1e(iTn),K2e(iTn),trange(iTn),'DIC'); %CO2aq,e at new T, from measured DIC
    RCint2=(Rnew*Rm2+ReG3)/(Rnew+1);
    [~,~,~,~,Roto3,~]=C13speciation(RCint,pHf(iTn),K1f(iTn),K2f(iTn),trange(iTn),'CO2aq');
    Doto3=(Roto3./stn-1).*1000;
    PcT{j,1}(k,17)=Doto3; %new otolith composition at varied T, constant d13met
    PcT{j,1}(k,18)=dm2; %new shared metabolic composition

    clearvars Tnew iTk2 Cr2 RCint Wnew ReG2 Rnew iTn Rint2 Doto2 Rm2 ReG3 RCint2 Roto3 Doto3 dm2;
    
    clearvars RoArag RoCO2aq RoCO2 dmod iTk tRef tRmf ttmpRe ttmpRo;
    end %length of data vector in each dataset
    
    %Generate model for mean dataset characteristics
    %To avoid weirdness with multiple datasets in each type (2 wild, 2
    %reared), just pick one temperature & type group for each model curve
    if j==1 %reared, 10C
        iw=[27:1:45,68:1:71]; %reared 10C
    else %wild, 4.5C
        %iw=[16:1:29]; %wild 4.5C
        %iw=[30:1:39]; %wild 2C from Jamieson only (not Gao, with different endmembers)
        iw=[40:1:49]; %wild 2C from Gao only (not Jamieson, with different endmembers)
    end
    
    tRef=(mean(d_dic(iw))./1000+1).*stn;tRmf=(mean(d_om(iw))./1000+1).*stn; %datum specific carbon endmembers
    tmean=mean(t(iw));
    [~,iTk]=min(abs(trange-tmean)); %approximate mean temperature, allows use of existing coefficients on trange  
    [~,ReG,~,~,~,~]=C13speciation(tRef,pHe(iTk),K1e(iTk),K2e(iTk),trange(iTk),'DIC'); %CO2aq,e   
    ttmpRe=ReG;clearvars ReG;
    %Temperature only model
    RoCO2=(tRmf.*Vh.*Pcvec./pCO2ref.*qd./qsT(iTk).*exp(-Eo(iTk).*TT(iTk))+ttmpRe)./...
        (Vh.*Pcvec./pCO2ref.*qd./qsT(iTk).*exp(-Eo(iTk).*TT(iTk))+1);
    CCvec=Vh.*Pcvec./pCO2ref.*qd./qsT(iTk).*exp(-Eo(iTk).*TT(iTk)); %vector of Cm/Ce that goes with this model
    for i=1:length(Pcvec)
    [~,~,~,~,R13Arag,~]=C13speciation(RoCO2(i),pHf(iTk),K1f(iTk),K2f(iTk),trange(iTk),spbC);
    RoArag(1,i)=R13Arag;
    end
    clearvars R13Arag RoCO2;   
    PcT{j,2}(:,1)=RoArag;
    PcT{j,2}(:,2)=Pcvec;
    PcT{j,2}(:,3)=tmean;
    PcT{j,2}(:,4)=trange(iTk);
    PcT{j,2}(:,5)=CCvec;
    
    clearvars RoArag RoCO2aq RoCO2 dmod iTk tRef tRmf ttmpRe ttmpRo CCvec;
    clearvars  t b d_dic d_dom d_meas;
end
%% Plotting for activity level
%2C to 14C
%2,4,4.5,7,10,14
clrsp=[0 0 0.9;...
    0.3 0.3 1;...
    0.5 0.5 1;...
    1 0.5 0;...
    1 0.3 0.3;...
    0.9 0 0];

% NOTE: There are 3 wild cod data that give large positive or negative
% values depending on the precise value of e_om (offset between blood
% dic and diet/tissue/otolith organic matter), because the denominator of
% the isotope ratio term (R above) is very close to zero:
% index values 26 and 28 have deloto-delmet=0.02 (same as best typically 
%obtained precision of individual carbonate measurements), index value 30 
% has deloto-delmet=-0.35. In all three cases this could be a problem with 
%the mean e_om chosen or the measurement of organic matter composition in 
%the underlying data, which is more variable and less certain within an 
%animals tissues than the precision of an individual carbonate sample
%(former is order of 0.1-0.4 replicate precision from same organic matter 
%source depending on quantity of material examined, similar to replicate 
%precision across different animals in same population). %For these very 
%small differences (<abs(~0.7)), accurate PhiCrit cannot be determined. See
%uncertainty analysis at end of script.

%Plot against Cr/Ce 
subplot(224); %pCO2m/pCO2e
%first plot model lines
plot(PcT{1,2}(:,5),(PcT{1,2}(:,1)./stn-1).*1000,'-','Color',clrsp(5,:),'LineWidth',2);hold on; %10C
plot(PcT{2,2}(:,5),(PcT{2,2}(:,1)./stn-1).*1000,'-','Color',clrsp(1,:),'LineWidth',2); %2C

%plot data
for jj=1:2
    xx=PcT{jj,1}(:,5); %pCO2m/pCO2e from sample
    yy=PcT{jj,1}(:,2); %Measured Ro_arag
    t=PcT{jj,1}(:,4); 
    itx=PcT{jj,1}(:,3); 

    if jj==1 %reared
        for i=1:length(t)
            if t(i)==2 clr=clrsp(1,:);
            elseif t(i)==4 clr=clrsp(2,:);
            elseif t(i)==4.5 clr=clrsp(3,:);
            elseif t(i)==7 clr=clrsp(4,:);
            elseif t(i)==10 clr=clrsp(5,:);
            elseif t(i)==14 clr=clrsp(6,:);
            end
            if xx(i)>=0
        plot((xx(i)),yy(i),'o','MarkerSize',4.5,'MarkerEdgeColor',clr,'MarkerFaceColor',clr,'LineWidth',0.1);hold on;
            end
        end
    else %wild
        for i=1:length(t)
            if t(i)==2 clr=clrsp(1,:);
            elseif t(i)==4 clr=clrsp(2,:);
            elseif t(i)==4.5 clr=clrsp(3,:);
            elseif t(i)==7 clr=clrsp(4,:);
            elseif t(i)==10 clr=clrsp(5,:);
            elseif t(i)==14 clr=clrsp(6,:);
            end
            if xx(i)>=0
        plot((xx(i)),yy(i),'o','MarkerSize',4.5,'MarkerEdgeColor',clr,'MarkerFaceColor','w','LineWidth',1.5);hold on;
            end
        end
    end
end
xlabel('P^{C}_{met} / P^{C}_{w}');ylabel('\delta^{13}C_{o}');
xlim([0 10]);
ylim([-6 0]);
set(gca,'Layer','top');
%% Histogram of PhiCrit
%Plot against Cr/Ce instead
subplot(223); %pCO2m/pCO2e
xx1=PcT{1,1}(:,1); %Estimated PhiCrit from reared
xx2=PcT{2,1}(:,1); %Estimated PhiCrit from wild
iflier=abs(xx2)>10; %see note below
xx2(iflier)=[];
%raw mean and std, mean and std of lognornal distribution
[mean(xx1),std(xx1),mean(xx2),std(xx2);...
    mean(log(xx1)),std(log(xx1)),mean(log(xx2)),std(log(xx2))] 
%mean and variance (not std) of this this lognormal distribution in linear
%space
[m1,v1]=lognstat(mean(log(xx1)),std(log(xx1)));
[m2,v2]=lognstat(mean(log(xx2)),std(log(xx2)));
[m1,v1,m2,v2] %close but not identical to normal stats fit to data
hhh=[mean(log(xx1)),std(log(xx1)),mean(log(xx2)),std(log(xx2))]
[hks,pks,tkk]=kstest2(xx1,xx2)
[length(~isnan(xx1)),length(~isnan(xx2))]

%Fit distribution with pdf kernel
kn1=0.1;
kn2=kn1; %kernal for nonparametric pdf, bandwidth
ptmp1=fitdist(xx1,'Kernel','Kernel','normal','Support','positive','Width',kn1);
ptmp2=fitdist(xx2,'Kernel','Kernel','normal','Support','positive','Width',kn2);
dx=0.5;edgetmp=[0:dx:8];centtmp=[(0+dx/2):dx:(8-dx/2)];
ctmp1=cdf(ptmp1,edgetmp);ctmp2=cdf(ptmp2,edgetmp);
%Below is same as probability=pdf*dx, but flexible for changing increment of dx
%considered if variable bin sizes.
clearvars Ptmp1 Ptmp2;
for i=1:(length(edgetmp)-1) %calculate probabilities associated with every increment
    Ptmp1(i,1)=ctmp1(i+1)-ctmp1(i);
    Ptmp2(i,1)=ctmp2(i+1)-ctmp2(i);
end

binedge=edgetmp;
[n1,~]=histcounts(xx1,binedge);
[n2,~]=histcounts(xx2,binedge);
s1=n1./sum(n1);s2=n2./sum(n2);
for bi=1:(length(binedge)-1)
    xc=[binedge(bi)+0.01,binedge(bi+1)-0.01,binedge(bi+1)-0.01,binedge(bi)+0.01];
    yc1=[0,0,s1(bi),s1(bi)];yc2=[0,0,s2(bi),s2(bi)];
    if s1(bi)<s2(bi)&s1(bi)>0&s2(bi)>0
        fill(xc,yc2,'w','EdgeColor','r');hold on;fill(xc,yc1,'w','EdgeColor','k');
    elseif s1(bi)>s2(bi)&s1(bi)>0&s2(bi)>0
        fill(xc,yc1,'w','EdgeColor','k');hold on;fill(xc,yc2,'w','EdgeColor','r');
    elseif s1(bi)>0&s2(bi)==0
        fill(xc,yc1,'w','EdgeColor','k');hold on;
    elseif s2(bi)>0&s1(bi)==0
        fill(xc,yc2,'w','EdgeColor','r');hold on;
    end
end

%Plot distributions
plot(centtmp,Ptmp1,'-','LineWidth',1.5,'Color','k');hold on;
plot(centtmp,Ptmp2,'-','LineWidth',1.5,'Color','r');
plot([0 8],[0 0],'-k','LineWidth',1);
set(gca,'TickLength',[0.015 0.015],'YAxisLocation','left');
xlabel('\Phi_{crit}');ylabel('{\itf}_{observed}');
xlim([0 6]);%ylim([0 0.17]);
xticks([0:1:6]);set(gca,'Layer','top');

%% Sensitivity of diagnostic ratio PCmet/PCw to uncertainty in measurements
%The ratio PCmet/PCw increases to infinity as the difference between
%del13Cint and del13Cmet decreases to zero. Thus analytical errors or poor
%estimates of either can lead to very large (+ or -) estimates of
%PCmet/PCw that are not meaningful. This section evaluates the sensitivty
%of the solution to variance in those two parameters in order to evaluate
%whether the empirical threshold identified for del13Cint-del13Cmet<0.5 to 
%0.75 is justified.

%Consider just replicate uncertainty (e.g., not systematic uncertainty in
%estimating relationships of diet, tissue, and metabolic carbon, chemical
%speciation and fractionation, poorly constrained environmental
%conditions) as minimum bound on random errors. 

%Metabolic carbon
%Gao study has tightly controlled diet and experimental conditions.
sdm(1,1)=std(C_IS{5,1}(1:63,5)); %all animals with measurements
sdm(2,1)=std(C_IS{5,1}(1:11,5));%different temperature groups
sdm(3,1)=std(C_IS{5,1}(12:26,5));
sdm(4,1)=std(C_IS{5,1}(27:45,5));
sdm(5,1)=std(C_IS{5,1}(46:63,5));
%Jamieson study has wild animals with similar measurements
sdm(1,2)=std(C_IS{6,1}(:,5)); %all animals with measurements
sdm(2,2)=std(C_IS{6,1}(1:15,5));%different temperature groups
sdm(3,2)=std(C_IS{6,1}(16:29,5));
sdm(4,2)=std(C_IS{6,1}(30:end,5));
sdm(5,2)=nan;
sda(1,:)=nanmean(sdm(2:end,:),1);
%Environmental replicates not more variable than reared, and both similar to the
%uncertainty in analytical replicates for the 13C of organic matter, and
%generally similar to or less than biological replicates typically reported
%for other fish species (order of 0.2 to 2 per mil).

%Otolith carbon (hence internal carbon)
%Evaluating Gao otolith compositions again provides biological replicate 
% uncertainty under highly controlled conditions.
sdo(1,1)=std(C_IS{5,1}(1:63,1)); %all animals with measurements
sdo(2,1)=std(C_IS{5,1}(1:11,1));%different temperature groups
sdo(3,1)=std(C_IS{5,1}(12:26,1));
sdo(4,1)=std(C_IS{5,1}(27:45,1));
sdo(5,1)=std(C_IS{5,1}(46:63,1));
%Jamieson study has wild animals with similar measurements, though range of
%activities means these might not be less similar.
sdo(1,2)=std(C_IS{6,1}(:,1)); %all animals with measurements
sdo(2,2)=std(C_IS{6,1}(1:15,1));%different temperature groups
sdo(3,2)=std(C_IS{6,1}(16:29,1));
sdo(4,2)=std(C_IS{6,1}(30:end,1));
sdo(5,2)=nan;
sda(2,:)=nanmean(sdo(2:end,:),1);
%The wild animals are consistently more variable than the reared, as
%expected given the wide range of ecological activity inferred for the wild
%animals. The reared animals are not far from the expected analytical
%uncertainty of ~0.2 per mil for carbonates.

%Seawater DIC (hence CO2)
%Gao, highly controlled, though with net pens environmental variability
%over time likely larger than reflected in reported values. Note that DIC
%not separately measured for each sample, generally only 2 samples per
%temperature group, may need to pool all even though d13C(DIC) varies with 
% temperature.
sdw(1,1)=std(C_IS{5,1}([1,7,12,17,27,36,46,59,78],4)); %all animals with measurements
sdw(2,1)=std(C_IS{5,1}([1,7,78],4));
sdw(3,1)=std(C_IS{5,1}([12,17],4));
sdw(4,1)=std(C_IS{5,1}([27,36],4));
sdw(5,1)=std(C_IS{5,1}([46,59],4));
%%No info on variance from Jamison
sdw(1:5,2)=nan;
sda(3,:)=nanmean(sdw(2:end,:),1);
%This seawater variance is lower than the typical replicate
%uncertainty for DIC samples, ~0.1-0.2 per mil, but the way it is reported
%(a single measurement for each group of concurrent samples) suggests
%replicate uncertainty was not fully captured.

%What is expected error in quadrature for ratio? (substituting 0.2 per mil for not
%measured wild DIC)
sdwtmp=sdw;sdwtmp(:,2)=0.2;
sn=sqrt(sdo.^2+sdwtmp.^2);
sd=sqrt(sdo.^2+sdm.^2); %0.4 to 1 per mil
sr=sqrt(((sn./(-10+15)).^2)+((sn./(-15+20)).^2));
%this is relative to value of ratio, i.e., 
% d(r)/r=sqrt([(ssn)/(dw-di)]^2+[(ssd)/(di-dm)]^2), so increases as the
% ratio does. The above example is just for hypothethical measurements 
%where the ratio of PCmet/PCw=1, giving uncertainty on the order of 10%.
% This will be easier to evaluate systematically using a Monte Carlo
% approach.

%Considering just the denominator term which dominates the uncertainty, the
%standard deviation is 0.4 to 1 per mil; considering just the reared
%experiments where activity variation is small and otolith composition is
%closer to a biological replicate for the same conditions gives sd=0.70 per
%mil. Thus differences<0.7 per mil are often within one standard deviation 
%of 0, suggesting this is a reasonable cutoff to avoid estimates of 
% PCmet/PCw which cannot be accurately constrained from measurements.

%Considering how the error changes as a function of the ratio for an
%idealized vector of isotopic values and replicate uncertainties from the
%reared dataset:
dw=-10;dm=-20;di=[(-10:-0.1:-20)];
RMm=sda(1,1).*randn(100000,length(di)); %Monte Carlo estimates around mean
RMi=sda(2,1).*randn(100000,length(di)); %using sd values across reared groups, not wild
RMw=sda(3,1).*randn(100000,length(di));
Rn=((dw+RMw)-(di+RMi))./((di+RMi)-(dm+RMm)); %Ratio with random errors
Ra=(dw-di)./(di-dm); %prescribed ratio
mm=prctile((Rn),50,1); %mean is biased by a few outliers, so use 50th percentile
ml=prctile((Rn),25,1);
mh=prctile((Rn),75,1);

figure;
ax=gca;
yl=[0 20];
p1=plot(di-dm,ml,'-','Color',[1 0.6 0.6]);hold on; %10th and 90th percentiles of errors
p2=plot(di-dm,mh,'-','Color',[1 0.6 0.6]);
p3=plot(di-dm,mm,'-r');
p4=plot(di-dm,Ra,'-k');
xlabel('\delta^{13}C_{int}-\delta^{13}C_{met}');
ylabel('P^{C}_{met}/P^{C}_{w}');
ylim(yl);

%What would pOw threshold scale be?
pO=(pCO2ref).*yl./(qd./qsT(iT));hold on;
yyaxis right;
p5=plot([0 10],[0.209 0.209],'-b');
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
ylim(pO);
ylabel('Hypoxia threshold (P^{O}_{w}, atm)');
%what is high end estimate of climatological pO2 (since supersaturation can
%be consistent in some regions)
%maxpo2=0.4486; %maximum climatological pO2 from WOA
maxpo2=0.2724; %99% instead
p6=plot([0 10],[maxpo2 maxpo2],'--b');
p1.ZData = ones(size(p1.XData));p2.ZData = ones(size(p2.XData));p3.ZData = ones(size(p3.XData));p4.ZData = ones(size(p4.XData));
p5.ZData = zeros(size(p5.XData));p6.ZData = zeros(size(p6.XData));
set(gca,'SortMethod','depth');
xticks([0:1:10]);

%This sensitivity analysis confirms that the ratio of PCmet/PCw becomes
%increasingly subject to large random errors as the difference approaches
%zero. The interquartile range of the expected errors covers most of the
%ocean pO2 range, such that by del13int-del13met<0.7 per mil, >25% of data 
% would be expected to exceed physical bounds (e.g., require hypoxia 
% threshold higher than could be supported by ocean oxygen; however a smaller 
% than average PCw would raise the POw associated with a given ratio). So
% the expected 1-sd analytical errors correspond to both increasingly
% unconstrained ratios and high likelihood of unphysical results. As the
% value gets even closer to zero, uncertianties on the three parameters
% contribute to a a systematic negative bias in the ratio because of
% spurious negative differences in the numerator or denominator when either
% measurement pair is within analytical uncertainties of zero difference.
 
% With the assumed endmembers above, this corresponds to a
%ratio of PCmet/PCw~13. Based on the database in Table S1, all but one
%experiment have a ratio of <14. One outlier, a pollock "exhaustively excercised"
%in a small tank had a ratio of ~35; this blood result could be real and  reflect
% induced anaerobic metabolism, or an systematic error from assuming too 
% low a value of PCw (and thus calculating too high a ratiovfrom PCint). 
% Regardless such a ratio could not be reliably measured using otoliths, 
% as the uncertainty bounds would be similar to the magnitude of the 
% estimate and encompass zero.
% 
%The above considerations suggest PCmet/PCw cannot be meaningfully constrained
% below del13int-del13met~0.7 per mil. Uncertainties become statistically
% limiting (e.g., 95% CIs encompass 0 to double the estimate) at even higher
% thresholds, so higher cutoffs would be more conservative with respect to 
% excluding potentially meaningless estimates, at the cost of also excluding 
%potentially valid results that could be grossly informative with appropriately large
%error bounds.
%
%In summary using a del13int-del13met cutoff of 0.7 per mil to screen out
%poorly constrained data and likely "fliers" is reasonable, and higher 
%thresholds could be justified. This sensitivity analysis also suggests 
%that uncertainty in estimated PCmet/PCw increases with the value of this
%ratio.
