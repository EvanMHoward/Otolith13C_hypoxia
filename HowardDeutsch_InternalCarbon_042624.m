%Author: Evan Howard (ehoward2@uw.edu)
%Original version: February 2022
%Current version: December 2022

%Description: This script evaluates distributions of internal carbon, and
%considers the relative roles of different isotopic endmembers on otolith
%compositions based on model sensitivities.

%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 

%% Prepare workspace
clear all;close all;clc;
addpath('.../Helper_Scripts_and_files/'); %Location of helper scripts

%Default parameters
%Van't Hoff equation constants
kb=8.617e-5;T0=273.15;Tref=15;
%endmembers
e_om=0; %approximate fractionation of blood DIC versus diet/tissue organic matter
        %if set to 0, d13met = d13om for plotting, but need to include
        %offset if converting real d13om data to d13met
del_dic=0; %d13C DIC 
del_om=-20+e_om; %d13C organic matter in diet or tissue
%environmental chemistry
pCO2_meas=370e-6; %approximate annual mean upper ocean pCO2
pO2_meas=0.209; %approximate annual mean upper ocean pO2
spbC='CO2aq'; %blood carbon species for mass balance
%animal
Pc=2.5; %default PhiCrit, active to resting oxygen supply to demand ratio
eps=-0.05; %default MI allometric sensitivity exponent
Eo=0.3; %default temperature sensitivity exponent
dEdT=0.01; %temperature dependence of Eo (typical intraspecies value is 0.01)
ET=15; %temperature at which standard Eo is observed

Vh=0.05; %efault hypoxia vulnerability (atm) at 15 degrees (C)
Bref=200; %Default reference mass (g)
%ranges of temperature and mass to consider          
trange=0:0.5:30; %Temperature range considered, leaving salinity fixed
brange=logspace(1,4,100); %Mass range considered, 10g to 10 kg
brangeN=brange./Bref; %Normalized mass compared to reference

%switches
swF=1; %use different blood chemistry than seawater chemistry
swStoi=1; %add factors to qd or qs 
        fd=0.5; %Vh_old*0.5, adjusts Vh from assumed pO2_blood=0 to average
        %pO2_blood=0.5*pO2_environment
        fs=1; %factor of CO2 ventilation related to metabolism
swEoT=1; %Include approximated temperature sensitivity of Eo based on 
        %transition from ventilation limited to diffusion limited supply
        %(only matters if dEdT nonzero)
swContCe=1; %Replace line plots comparing Ce, Cm, Phi with second contour plot

if swEoT==1
    Eo=Eo+dEdT*(Tref-ET); % Eo value at Tref
end

%% Aqueous chemistry and switches
%Arrhenius temperature term
TT=(1/kb).*(1./(trange+T0)-1./(Tref+T0));

%Environmental gas partial pressures
pCO2ref=pCO2_meas;pO2ref=pO2_meas;

%Assumed pH and salinity/ionic strength of environment and fish fluids
pH_oc=8.07;TpHref_oc=10; %Reference pH (total scale) and temperature, based on GLODAP
pH_blood=7.85-0.15;TpHref_blood=TpHref_oc; %pH_NBSscale of fish blood as 
%low as 7.3, but generally around 7.7-8.0; -0.15 approximates scale difference
S_oc=34;S_f=9; %approximate ocean salinity, as well as S with equivalent 
%ionic stregth of fish plasma or endolymph
if swF~=1
    S_f=S_oc;pH_blood=pH_oc; %overwrite fish blood with seawater chemistry
end

%Temperature sensitivity of pH
pHe=pH_T(trange,S_oc,pH_oc,TpHref_oc);pHe=pHe';
pHf=pH_T(trange,S_f,pH_blood,TpHref_blood);pHf=pHf';
%pHe(:)=pH_oc;pHf(:)=pH_blood;

%Carbon dioxide solubility
[K0e K1e K2e]=carbeq(trange,S_oc);K0e=K0e';K1e=K1e';K2e=K2e'; %In seawater
[K0f K1f K2f]=carbeq(trange,S_f);K0f=K0f';K1f=K1f';K2f=K2f'; %In fish bodily fluids

%Oxygen solubility
Khe=O2sol(S_oc,trange)./(pO2ref).*1e-6;Khe=Khe'; %mol kg-1 atm-1 of O2
Khf=O2sol(S_f,trange)./(pO2ref).*1e-6;Khf=Khf'; 

%Diffusivities in water (all in 10-4 m2 s-1)
[Dco2e,Dhco3e,Dco3e]=co2_diffusion(S_oc,trange);Dco2e=Dco2e';Dhco3e=Dhco3e';Dco3e=Dco3e';
[Dco2f,Dhco3f,Dco3f]=co2_diffusion(S_f,trange);Dco2f=Dco2f';Dhco3f=Dhco3f';Dco3f=Dco3f';
[~,~,~,~,Do2e,~,~]=gas_diffusion(S_oc,trange);Do2e=Do2e';
[~,~,~,~,Do2f,~,~]=gas_diffusion(S_f,trange); Do2f=Do2f';

%Isotopic endmembers
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
qC=K0e.*Dco2e;
if swStoi==1
qs=qs.*fs;qd=qd.*fd;
end

%Temperature dependence of Eo
if swEoT==1
Eo=Eo+dEdT.*(trange-ET);
end

%% PhiCrit*pCO2crit/pCO2w vs Rm
%What temperature is this for?
load TnormCodOto_122922; %PcT
Tnew=PcT{1,1}(1,11);  %Tref or, use Tnew=PcT{1,1}(1,11) if using temperature corrected, filled data in scatter overlay
[~,iT2]=min(abs(trange-Tnew));

%figure graphical properties
figure('Renderer', 'painters', 'Position', [100 100 400 600]);

clearvars RoT;
qsT=qs(iT2); %fixed T and reference temperature
%Compare d13Co to Rm for various values of PhiCrit*pCO2crit/pCO2w
pCPcvec=[0.0:0.01:0.5].*qd./qsT./pCO2ref; 
%Based on compilation in Deutsch et al. 2020 (Nature), Vh varies from 0.003
%to 0.161 atm O2, PhiCrit from 1.1 to 7.4, and Vh*PhiCrit from 0.015 to 0.40.
dmfvec=[-25:0.5:-15]; %typical range for marine foodwebs, from 
%trophic levels 2-5 (e.g., planktiovres to top predators).

Rmfvec=((dmfvec+e_om)./1000+1).*stn; %variable diet endmember, including 
%any enrichment in blood DIC versus tissues.
for i=1:length(dmfvec)
RoCO2=(Rmfvec(i).*pCPcvec+Reref)./(pCPcvec+1);
[~,~,~,~,R13Arag,~]=C13speciation(RoCO2,pHf(iT2),K1f(iT2),K2f(iT2),trange(iT2),spbC);
RoT(i,:)=R13Arag;clearvars R13Arag RoCO2;
end
h222=subplot(211); %PhiCrit/pCO2 vs Rm
%Changed here to output directly as PhiCrit (at fixed value of Vh) instead of Pmet/Pw
Vhtmp=0.0621;Pctmp=pCPcvec./(qd./qsT./pCO2ref.*Vhtmp);%mode of all species Vh ~0.05, Atlantic cod 0.0621
[X,Y]=meshgrid(dmfvec,Pctmp);
aa=(RoT./stn-1).*1000;
contourf(Y,X,aa',[-12:1:2],'edgecolor','none');
%%%
colormap(h222,flipud(parula(14)));
cq2=[-12 2];clrc=cq2;c=colorbar('XTickLabel',{'','','',''});caxis(clrc);
c.XTick=[(cq2(1)+1):2:cq2(2)];c.TickLabels=[(cq2(1)+1):2:cq2(2)];
hold on;
ylabel('\delta^{13}C_{om}');
xlabel('\Phi_{c}');%xlabel('P^{C}_{met} / P^{C}_{w}');
xlim([0 8]);

%% Overylay measured cod d13Coto, normalized to constant temperature
for j=1:2 %reared, wild
    if j==1 cl='k';
        ibad=zeros(length(PcT{j,1}(:,1)),1);
    else cl='r';
        ibad=zeros(length(PcT{j,1}(:,1)),1);
        ibad(PcT{j,1}(:,6)<0.7)=1;
%     See explanation in HowardDeutsch_AtlanticCodCase_042624, 
%     can't determine accurately when internal and metabolic compositions within uncertainties.
    end
    cf=PcT{j,1}(~ibad,10); %Fill by estimated Coto if temperature normalized (10, or 2 if as measured)
    cf=round(cf*2,0)/2; %change to same discrete colorbar, rounded to nearest 0.5
    %use above as fill if there is a sensible way to temperature normalize,
    %but without knowing local water pCO2 this will lead to systematic biases 
    %between each group (particularly across the wild datasets). One
    %could estimate pCO2 based on local conditions plus solubility as in
    %the global analysis, but this is not needed to make the
    %point about modes of variations.
    hold on;scatter(PcT{j,1}(~ibad,1),PcT{j,1}(~ibad,7),45,cf,'filled','MarkerEdgeColor',cl,'LineWidth',1); 
    %NOTE: Cod includes d13met=d13om+1 per mil, here compared against model
    %that has d13met=d13om (so y axes do align as d13met)
end

%% PhiCrit*pCO2crit/pCO2w vs Re (as contours over Cm plot)
hold on;

if swContCe==1
clearvars RoT;
qsT=qs(iT2); %fixed T and reference temperature
defvec=[-15:0.5:-5]; %Replace with d13CO2aq
dm=PcT{1,1}(1,18);rm=(dm/1000+1)*stn;

Refvec=((defvec)./1000+1).*stn; %variable seawater endmember
for i=1:length(defvec)
[~,tmpRevec(i),~,~,~,~]=C13speciation(Refvec(i),pHe(iT2),K1e(iT2),K2e(iT2),trange(iT2),'CO2aq'); %CO2aq,e 
end
%Rmf is ratio for d13met, here assume d13met=d13om=-20
for i=1:length(defvec)
RoCO2=(rm.*pCPcvec+tmpRevec(i))./(pCPcvec+1);
[~,~,~,~,R13Arag,~]=C13speciation(RoCO2,pHf(iT2),K1f(iT2),K2f(iT2),trange(iT2),spbC);
RoT(i,:)=R13Arag;
end
clearvars R13Arag RoCO2 tmpRevec;

aa=(RoT./stn-1).*1000;

ax2=subplot(212);
[X,Y]=meshgrid(defvec,Pctmp);
contourf(Y,X,aa',[-6:1:2],'edgecolor','none');
colormap(ax2,flipud(parula(8)));
cq2=[-6 2];clrc=cq2;c=colorbar('XTickLabel',{'','','',''});caxis(clrc);
c.XTick=[(cq2(1)+1):2:cq2(2)];c.TickLabels=[(cq2(1)+1):2:cq2(2)];
ylabel('\delta^{13}C_{w}');
xlabel('\Phi_{c}');%xlabel('P^{C}_{met} / P^{C}_{w}');
xlim([0 8]);

%Overlay otolith data normalized to constant T and d13met, but measured
%d13w
for j=1:2 %reared, wild
    if j==1 cl='k';
        ibad=zeros(length(PcT{j,1}(:,1)),1);
    else cl='r';
        ibad=zeros(length(PcT{j,1}(:,1)),1);
        ibad(PcT{j,1}(:,6)<0.7)=1;
    end
    cf=PcT{j,1}(~ibad,17); %Fill by estimated Coto if temperature normalized (17, or 2 if as measured)
    cf=round(cf*2,0)/2; %change to same discrete colorbar, rounded to nearest 0.5
    hold on;scatter(PcT{j,1}(~ibad,1),PcT{j,1}(~ibad,8),45,cf,'filled','MarkerEdgeColor',cl,'LineWidth',1); 
end
end

%% O2 and CO2 based PhiCrit
figure('Renderer', 'painters', 'Position', [100 100 400 600]);
h2=subplot(211);
swPDF=1;

qsT=qs(iT2); %fixed T and reference temperature
qCT=qC(iT2);

%From respirometry
load Deutsch2020Phist.mat %2020 Nature paper traits for animals with pO2crit &/or PhiCrit
%column 1 Vh (atm O2), column 2 PhiCrit (90%), column 3 PhiCrit (F1)
%The two estimates of PhiCrit are usually similar enough that an average is
%OK for plotting a qualitative histogram here.
aa=Deutsch2020Phist(:,1).*0.5./qsT./pCO2ref; %Factor of 0.5 depends on whether 
%trying to match arterial, venous, or whole-body fluids (0.5 gives whole
%body average)
bb=nanmean(Deutsch2020Phist(:,2:3),2);
cc=aa.*bb;
%Fit PhiCrit to start
kn1=(1/5); %kernal for nonparametric pdf, for active
dx=0.5;edgetmp=[0:dx:20];centtmp=[(edgetmp(1)+dx/2):dx:(edgetmp(end)-dx/2)];
ptmpPhiCmeas=fitdist(bb,'Kernel','Kernel','normal','Support','positive','Width',kn1);
ctmpPhiCmeas=pdf(ptmpPhiCmeas,centtmp);
PtmpPhiCmeas=ctmpPhiCmeas.*(dx); %calculate probability assuming constant dx

%From blood pCO2 measurements
load BloodPint.mat;
A=Pint{:,[8,9]};AT=Pint{:,7};
Braw=(A(:,1)-A(:,2))./A(:,2); %as Pmet/Pw, normalized to different Pw values already
%find closest value of qs for experiments, use to normalize to constant
%temperature if needed
for i=1:length(AT)
[~,iTexp(i,1)]=min(abs(AT(i)-trange));
end
%Active and resting experiments are already well paired in T, so no need to
%scale for temperature differences.
B=Braw;

%All active
ih=zeros(size(Pint,1),1);ih([6;7;18;20;22;24;26;28;30;32;34;36;38;40;42;...
    44;46;48;50;52])=1;ih=logical(ih);
cc2=B(ih);aa2=B(~ih);
%Paired resting and active only
ir=[8,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51];
ia=[7,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52];%don't pair 5,6, Pw too different
bb2=B(ia')./B(ir');
%Fit distribution
ptmpPhiCmeas2=fitdist(bb2,'Kernel','Kernel','normal','Support','positive','Width',kn1);
ctmpPhiCmeas2=pdf(ptmpPhiCmeas2,centtmp);
PtmpPhiCmeas2=ctmpPhiCmeas2.*(dx); %calculate probability assuming constant dx

subplot(211);
binedge=[0:dx:7];
[n1,~]=histcounts(bb,binedge);
[n2,~]=histcounts(bb2,binedge);
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
plot(centtmp,PtmpPhiCmeas,'-','LineWidth',2,'Color','k');
plot([0 8],[0 0],'-k','LineWidth',1.5);
xlim([0 8]);ylim([0 0.35]);

plot(centtmp,PtmpPhiCmeas,'-','LineWidth',2,'Color','k');
plot(centtmp,PtmpPhiCmeas2,'-','LineWidth',2,'Color','r');
ylabel('\Phi_{c}');
xlim([0 8]);ylim([0 0.35]);
yticks([0:0.05:0.35]);

%Distribution tests
%From different distribution?
[hhks,ppks,kk]=kstest2(bb,bb2) %KS test, hh=0 if can't reject null hypothesis of same distribution
[hhtt,pptt,~,tt]=ttest2(log(bb),log(bb2)) %t-test on log-transformed, hh=0 if can't reject null hypothesis of same mean
length(~isnan(bb)),length(~isnan(bb2))

%NOTE: must print as svg to keep vector format and transparency. Otherwise
%'eps' printed is not a real eps but a bitmap (if transparent) or if a true
%eps has no transparency.

%% Resting Pmet comparison
swPmet=1;
    Eotmp=0.383; %mean Eo across database of Deutsch et al. 2020
if swPmet~=1
aa=Deutsch2020Phist(:,1).*0.5./qsT/pCO2ref;%.*1e6; %Pmet from O2, in uatm, all at 15C
A=Pint{:,[8,9]};AT=Pint{:,7};
%find closest value of qs for experiments, use to normalize to same T
for i=1:length(AT)
[~,iTexp(i,1)]=min(abs(AT(i)-trange));
end
Braw=(A(:,1)-A(:,2))./A(:,2); 
B=Braw.*qs(iTexp)./qsT.*exp((-Eotmp/kb).*(1./(Tref+273.15)-1./(trange(iTexp)+273.15)));
ivh=ih;%ivh([2;19;31;37;45;49],1)=1; %optionally also remove venous blood samples along with active
aa2=B(~ivh);%Pmet from CO2, in uatm
aa3=A(~ivh,1);
kn1=(1/5); %kernal for nonparametric pdf, fraction of range, for active
kn3=(1/5);
dx2=0.5;edgetmp=[-dx2:dx2:20];centtmp=[(edgetmp(1)+dx2/2):dx2:(edgetmp(end)-dx2/2)];
dx3=500;edgetmp3=[-dx3:dx3:15000];centtmp3=[(edgetmp3(1)+dx3/2):dx3:(edgetmp3(end)-dx3/2)];

ptmp1=fitdist(aa,'Kernel','Kernel','normal','Support','positive','Width',kn1);
ptmp2=fitdist(aa2,'Kernel','Kernel','normal','Support','positive','Width',kn1);
ptmp3=fitdist(aa3,'Kernel','Kernel','normal','Support','positive','Width',kn3);
ctmp1=pdf(ptmp1,centtmp);ctmp2=pdf(ptmp2,centtmp);ctmp3=pdf(ptmp3,centtmp3);
clearvars Ptmp1 Ptmp2 Ptmp3;
Ptmp1=ctmp1.*dx2;Ptmp2=ctmp2.*dx2;Ptmp3=ctmp3.*dx3; %probabilities

%figure graphical properties
figure('Renderer', 'painters', 'Position', [100 100 600 300]);
subplot(211);
patch([300 2000 2000 300],[0 0 0.2 0.2],'blue','FaceAlpha',0.3,'EdgeColor','none'); hold on;
plot((centtmp3),Ptmp3,'-','LineWidth',2,'Color',[0.85 0 0]);
plot(([0 10000]),[0 0],'-k','LineWidth',1.5);
ylim([0 0.2]);xlim([0 7000]);
subplot(212);
plot((centtmp),Ptmp1,'-','LineWidth',2,'Color','k'); hold on;
plot((centtmp),Ptmp2,'-','LineWidth',2,'Color',[0.85 0 0]); 
plot(([0 20]),[0 0],'-k','LineWidth',1.5);
plot(([1 1]),[0 1],'-b','LineWidth',1.5);
ylim([0 0.2]);xlim([0 7]);
ylabel('{\itf} observed');xlabel('P^{C}_{met} / P^{C}_{w}');
set(gca,'Box','off');
[hPPks,pPPks,tPPkk]=kstest2(aa,aa2)
end

if swPmet==1
%PCmet normalized via both Q(T) and Pw, calculate associated PCint
aa=Deutsch2020Phist(:,1).*0.5./qsT.*1e6+(pCO2ref.*1e6); %Pint from O2 + Pw, in uatm, all at 15C
aa3=aa-(pCO2ref.*1e6); %Pmet
A=Pint{:,[8,9]};AT=Pint{:,7};
%find closest value of qs for experiments, use to normalize to same T
for i=1:length(AT)
[~,iTexp(i,1)]=min(abs(AT(i)-trange));
end
Braw=(A(:,1)-A(:,2)); %Pint-Pw
%normalize PCmet to 15C and Pw=pCO2ref=370 uatm, and calc Pint=Pmet+Pw
%Same steps as line 483, but scale experimental Pmet/Pw to a new reference
%Pw for equivalence to oxygen estimates (could alternatively rescale oxygen
%estimates to experimental T and O2, but then a range of conditions would
%be plotted rather than a distribution at fixed conditions).
B=Braw.*qs(iTexp)./qsT.*...
    exp((-Eotmp/kb).*(1./(Tref+273.15)-1./(trange(iTexp)+273.15))).*...
    (pCO2ref.*1e6)./A(:,2)+(pCO2ref.*1e6); 
 
%Note: This normalization scheme takes advantage of observing that species'
% Pmet/Pw is relatively stable and agrees well when derived from either O2 or CO2.

ivh=ih;%ivh([2;19;31;37;45;49],1)=1; %optionally also remove venous blood samples along with active
aa2=B(~ivh);%Pint from CO2, in uatm
aa4=aa2-(pCO2ref.*1e6); %Pmet
kn1=(1/5); %kernal for nonparametric pdf, fraction of range, for active
dx3=500;edgetmp3=[-dx3:dx3:15000];centtmp3=[(edgetmp3(1)+dx3/2):dx3:(edgetmp3(end)-dx3/2)];
kn3=(1/5);
ptmp1=fitdist(aa,'Kernel','Kernel','normal','Support','positive','Width',kn3);
ptmp2=fitdist(aa2,'Kernel','Kernel','normal','Support','positive','Width',kn3);
ptmp3=fitdist(aa3,'Kernel','Kernel','normal','Support','positive','Width',kn3);
ptmp4=fitdist(aa4,'Kernel','Kernel','normal','Support','positive','Width',kn3);
ctmp1=pdf(ptmp1,centtmp3);ctmp2=pdf(ptmp2,centtmp3);
ctmp3=pdf(ptmp3,centtmp3);ctmp4=pdf(ptmp4,centtmp3);
clearvars Ptmp1 Ptmp2 Ptmp3 Ptmp4;
Ptmp1=ctmp1.*dx3;Ptmp2=ctmp2.*dx3; %probabilities
Ptmp3=ctmp3.*dx3;Ptmp4=ctmp4.*dx3;

%figure graphical properties
figure('Renderer', 'painters', 'Position', [100 100 400 600]);

% histograms where distribution with lowest number in a bin appears on
%as the front image for that bin.
binedge=[0:500:6000];
subplot(211);
[n1,~]=histcounts(aa,binedge);
[n2,~]=histcounts(aa2,binedge);
s1=n1./sum(n1);s2=n2./sum(n2);
[n3,~]=histcounts(aa3,binedge);
[n4,~]=histcounts(aa4,binedge);
s3=n3./sum(n3);s4=n4./sum(n4);
for bi=1:(length(binedge)-1)
    xc=[binedge(bi)+5,binedge(bi+1)-5,binedge(bi+1)-5,binedge(bi)+5];
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
plot(([pCO2ref.*1e6 pCO2ref.*1e6]),[0 1],'-b','LineWidth',4);hold on;
plot((centtmp3),Ptmp1,'-','LineWidth',2.5,'Color','k'); hold on;
plot((centtmp3),Ptmp2,'-','LineWidth',2.5,'Color',[0.85 0 0]); 
plot(([0 20]),[0 0],'-k','LineWidth',1.5);
ylim([0 0.6]);xlim([0 4500]);%change plotted x limit here
ylabel('{\itf} observed');xlabel('P^{C}_{int}');
set(gca,'Box','off');
subplot(212);
for bi=1:(length(binedge)-1)
    xc=[binedge(bi)+5,binedge(bi+1)-5,binedge(bi+1)-5,binedge(bi)+5];
    yc1=[0,0,s3(bi),s3(bi)];yc2=[0,0,s4(bi),s4(bi)];
    if s3(bi)<s4(bi)&s3(bi)>0&s4(bi)>0
        fill(xc,yc2,'w','EdgeColor','r');hold on;fill(xc,yc1,'w','EdgeColor','k');
    elseif s3(bi)>s4(bi)&s3(bi)>0&s4(bi)>0
        fill(xc,yc1,'w','EdgeColor','k');hold on;fill(xc,yc2,'w','EdgeColor','r');
    elseif s3(bi)>0&s4(bi)==0
        fill(xc,yc1,'w','EdgeColor','k');hold on;
    elseif s4(bi)>0&s3(bi)==0
        fill(xc,yc2,'w','EdgeColor','r');hold on;
    end
end
plot(([pCO2ref.*1e6 pCO2ref.*1e6]),[0 1],'-b','LineWidth',4);hold on;
plot((centtmp3),Ptmp3,'-','LineWidth',2.5,'Color','k'); hold on;
plot((centtmp3),Ptmp4,'-','LineWidth',2.5,'Color',[0.85 0 0]); 
plot(([0 20]),[0 0],'-k','LineWidth',1.5);
ylim([0 0.55]);xlim([0 4500]);%change plotted x limit here
ylabel('{\itf} observed');xlabel('P^{C}_{met}');
set(gca,'Box','off');

%Check if no difference in distribution (null). KS test is more senstivie
%to central tendancy and less sensitive to tails.
[hPPks,pPPks,tPPkk]=kstest2(aa,aa2);H1(1,:)=[hPPks,pPPks,tPPkk];
[hPPdks,pPPdks,tPPdkk]=kstest2(Ptmp1,Ptmp2);H1(2,:)=[hPPks,pPPks,tPPkk];
[hPPks,pPPks,tPPkk]=kstest2(aa3,aa4);H1(4,:)=[hPPks,pPPks,tPPkk];

H1
size(aa),size(aa2)

end

