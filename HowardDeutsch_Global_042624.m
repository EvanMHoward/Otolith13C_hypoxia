%Author: Evan Howard (ehoward2@uw.edu)
%Original version: March 2022
%Current version: April 2024

%Description: This script generates distributions of otolith carbon
%isotopic ratios expected across distributions of hydrography, trophic
%level, and hypoxia traits. The outputs of this script (d3D.mat) can are
%further processed and plotted using the associated script
%'HowardDeutsch_GlobalCase_....m'. This field can be regenerated with other
%choices about how the distributions are used (tens of minutes), or simply 
%loaded from the archived fields. 

%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 

%% Prepare workspace
clear all; close all; clc;

%Current filepath
addpath('.../Helper_Scripts_and_files/'); %Location of helper scripts

%Default parameters
zmax=200; %integration depths for maps and figures
dEdT=0.0; %Good approximation, could do more complicated caluclation with 0.02 eV/C
%Trophic level relationships based on 15 species datasets, see SM.
t_om=-19; %OM-DIC offset at trophic level 4.1 (see HowardDeutsch_OMminDIC_corr_042624 in helper scripts)  
ts_om=2; %standard deviation of OM-DIC offset
TLd=4.1;
TLE=1.5; %approximate ecosystem average trophic enrichment per trophic level, see SM
e_om=0; %additional fractionation factor for metabolic endmember relative to tissue composition (if desired)
pO2atm=0.209; %atmospheric oxygen partial pressure
spbC='CO2aq'; %relevant blood carbon species from mass balance

%Switches
% selectively remove parameter variability to explore drivers (1=average, 0=variable)
domeanChem=0; % all water properties (except T, O2)
domeanTfrac=0; % T-dependent fractionation
domeanTsup=0; % T-dependent supply (Qs)
domeanTbio=0; % T-dependent Pcrit (same as setting Eo=0)
domeanTrait=0; % species traits
swF=1; %use different blood chemistry than seawater chemistry
swStoi=1; %add factors to qd or qs 
        fd=0.5; %Vh_old*0.5, adjusts Vh from assumed pO2_blood=0 to average
        %pO2_blood=0.5*pO2_environment
        fs=1; %factor of CO2 ventilation related to metabolism
swEoT=1; %Include approximated temperature sensitivity of Eo based on 
        %transition from ventilation limited to diffusion limited supply
        %(only affects result if dEdt nonzero.
swSave=0; %Save 3d output for plotting instead of regenerating each time
         %%%If 0, load existing processed fields.

%Van't Hoff equation constants
kb=8.617e-5;T0=273.15;Tref=15;

%% Load species datasets
%Respirometry based MI traits (Deutsch et al. 2020)
load Phy; %Has structures Phy (species-spepcific traits) and Gtr (fit distributions of traits)
clearvars Phy;

%Trophic level
load TrophicLevel; %FishBase trophic level database
FoodTroph=TrophicLevel{:,2};
i=isnan(FoodTroph);
FoodTroph(i)=TrophicLevel{1,3}(i);
%fill gaps with alternative trophic level metric, which is usually similar 
%at the individual level and indistinguishable statistically across the database.
DelOMDIC=t_om+TLE.*(FoodTroph-TLd);
%Fit distribution to trophic level
kn1=0.05; %kernal for nonparametric pdf, fraction of range
%e.g., if range is -23 to -18, then kn1=0.1 smooths over 0.5 per mil
pdfD=fitdist(DelOMDIC,'Kernel','Kernel','normal','Support','unbounded','Width',kn1);
DCent=[-23:0.25:-18]; %vector to evaluate pdf of OM-DIC trophic enrichment
pD=pdf(pdfD,DCent);
%pdf, NOT probability like histogram; convert to cdf and integrate over 
%specific range to evaluate associated probabilities if desired.
%figure;subplot(211);plot([-23:0.25:-18],pD,'-k');subplot(212);histogram(DelOMDIC,[-22.75:0.25:-18.25]);

% Otolith isotope data
load InterSpecies_b.mat

%% Ocean fields
vnames={'temp','salt','o2','po4','aou'};
load WOA; %World ocean atlas 2018
load GLO; %GLODAPv2 2016 climatology
load Goc; %WOA grids
%Climatological AOU
o2sat=O2sol(WOA.salt,WOA.temp);
WOA.aou=o2sat-WOA.o2;

% subset Ocean
WOAsub=subset_WOA(WOA,vnames,zmax);
Iz=WOAsub.Iz;
WOAsub.pco2=GLO.pco2(:,:,Iz);
WOAsub.dic=GLO.DIC(:,:,Iz);
WOAsub.alk=GLO.Alk(:,:,Iz);
WOAsub.pH=GLO.pH(:,:,Iz);
v3d=WOA.v3d(:,:,Iz); %ocean volume, m3, in each cell

%DIC/phosphate anomaly correlation
load mdl13DIC;b=mdl13DIC{1};
dic13=b(1)*WOA.temp+b(2)*WOA.salt+b(3)*WOA.po4+b(4)*WOA.aou+b(5);

aveT=nanmean(WOAsub.temp(:));
aveO=nanmean(WOAsub.po2(:));
aveD=nanmean(WOAsub.dic(:));
aveC=nanmean(WOAsub.pco2(:));
aveA=nanmean(WOAsub.alk(:));
aveS=nanmean(WOAsub.salt(:));
aveH=nanmean(WOAsub.pH(:));
aved=nanmean(dic13(:));

clearvars WOA GLO Vol;

%% Select traits and prepare for calculation
%Probability object vectors
dxw=0.2;bed=[(-25-dxw/2):dxw:(10+dxw/2)];bcent=[-25:dxw:10];
% Preallocate 3d arrays
wtPDFall=nan(size(WOAsub.pco2,1),size(WOAsub.pco2,2),size(WOAsub.pco2,3),length(bcent));
d13oto_3d_AEave=nan(size(WOAsub.pco2));
%d13oto_3d_AEstd=nan(size(WOAsub.pco2));
d13int_3d_AEave=nan(size(WOAsub.pco2));
%d13int_3d_AEstd=nan(size(WOAsub.pco2));
d13om_3d_DODave=nan(size(WOAsub.pco2));
%d13om_3d_DODstd=nan(size(WOAsub.pco2));
d13met_3d_DODave=nan(size(WOAsub.pco2));
%d13met_3d_DODstd=nan(size(WOAsub.pco2));
d13co2_3d=WOAsub.pco2*nan;
d13dic_3d=WOAsub.pco2*nan;
d13om_3d=WOAsub.pco2*nan;
d13met_3d=WOAsub.pco2*nan;
dic_3d=WOAsub.pco2*nan;
pco2_3d=WOAsub.pco2*nan;
PcritO2_3d=WOAsub.pco2*nan;
VhPhic_3d=WOAsub.pco2*nan;
Eo_3d=WOAsub.pco2*nan;
PcritC_3d=WOAsub.pco2*nan;

%Trait data
Trt.gA=Gtr.gA; % Values of Ac (=Ao/Phi_crit)
Trt.gE=Gtr.gE; % Values of Eo
Trt.AEdist=Gtr.AcEfit; % Trait pdf (Ac, Eo)
Trt.dEdT=dEdT; %dEdT as set in initial parameters
Trt.phimax=200; % max Phi (set >100 to be functionally unlimited)
Trt.tOMdist=pD; %Delta(OM-DIC) pdf
Trt.tOM=DCent; %Detla(OM-DIC) evaluation values

%% Calculate global otolith values
if swSave==1

% loop over spatial dimensions
for k=1:length(WOAsub.z)%1:length(WOAsub.z)
    [k length(WOAsub.z)]
    tic;
    for i=1:length(WOAsub.x)
        for j=1:length(WOAsub.y)
            
            % extract environment condition
            Env.T=nanmean(WOAsub.temp(i,j,k,:),4);
            Env.S=nanmean(WOAsub.salt(i,j,k,:),4);
                if isnan(Env.S); Env.S=34;end %approximate salinity for missing values
            Env.O=nanmean(WOAsub.po2(i,j,k,:),4); % atm
            Env.D=WOAsub.dic(i,j,k); 
            Env.C=WOAsub.pco2(i,j,k); 
            Env.A=WOAsub.alk(i,j,k); 
            Env.H=WOAsub.pH(i,j,k); 
            Env.d=dic13(i,j,k);
            Env.Tbio=Env.T; % T used in physiology
            Env.Tsup=Env.T; % T used in Qs (supply)
            Env.Tfrac=Env.T; % T used in C fractionation
            
            if domeanChem % overwrite with averages
                Env.D=aveD;Env.C=aveC;Env.A=aveA;Env.S=aveS;Env.H=aveH;Env.d=aved;
            end
            if domeanTfrac;Env.Tfrac=aveT;end
            if domeanTsup;Env.Tsup=aveT;end
            if domeanTbio;Env.Tbio=aveT;end
            
            if isfinite(Env.C + Env.H + Env.T + Env.S) % skip cells w/o data 
                jjj=1;
                if jjj==1
            %Arrhenius temperature term
            TT=(1/kb).*(1./(Env.Tbio+T0)-1./(Tref+T0));    
            %Environmental gas partial pressures
            pO2ref=pO2atm;pCO2ref=Env.C.*1e-6;
            %Fish fluid ionic composition
            pHf=(7.85-0.15);S_f=9;%Salinity equivalent of blood ionic strength for gas chemistry
            if swF~=1; pH_blood=Env.H;end %overwrite fish blood with seawater chemistry
            %No temperature-scaling of pH if assumed to be at measured ocean T
            %Carbon dioxide solubility
            [K0e K1e K2e]=carbeq(Env.Tsup,Env.S);[K0f K1f K2f]=carbeq(Env.Tsup,S_f);
            %Oxygen solubility
            Khe=O2sol(Env.S,Env.Tsup)./(pO2ref).*1e-6; Khf=O2sol(S_f,Env.Tsup)./(pO2ref).*1e-6;%mol kg-1 atm-1 of O2 
            %Diffusivities in water (all in 10-4 m2 s-1)
            [Dco2e,Dhco3e,Dco3e]=co2_diffusion(Env.S,Env.Tsup);[Dco2f,Dhco3f,Dco3f]=co2_diffusion(S_f,Env.Tsup);
            [~,~,~,~,Do2e,~,~]=gas_diffusion(Env.S,Env.Tsup);[~,~,~,~,Do2f,~,~]=gas_diffusion(S_f,Env.Tsup);
            %Supply stoichiometry
            qd=1;qs=K0e.*Dco2e./(Khe.*Do2e);
            if swStoi==1;qs=qs.*fs;qd=qd.*fd;end
            %Temperature dependence of Eo
            if swEoT==1;Eo=Trt.gE+Trt.dEdT.*(Env.Tbio-Tref);
            else Eo=Trt.gE;end
            %Isotopic endmembers
            stn=0.011180;
            def=Env.d;
            dmf=def+DCent; %this is 1xlength(DCent) vector
            Ref=((def)./1000+1).*stn;Rmf=((dmf)./1000+1).*stn;
            [~,Rew,~,~,~,~]=C13speciation(Ref,Env.H,K1e,K2e,Env.Tfrac,'DIC'); %CO2aq,seawater  
            
            %Biological rates (no biomass scaling)
            VhPhic=(1./Trt.gA); % Vh*PhiCrit, because gA=Ao/PhiCrit
            PcritO2active=VhPhic.*exp(-Eo.*TT); %Vh(T)*PhiCrit
            Phic=Env.O./PcritO2active; % ratio of local pO2/[active pO2 threshold],
            PcritO2active(Phic<1 | Phic>Trt.phimax)=nan; % WARNING - maximum Phi?
            PcritC=PcritO2active.*(qd./qs);
            
            %Internal fluid and otolith composition
            Rint=nan(size(PcritC,1),size(PcritC,2),length(Rmf));
            d13_int=nan(size(PcritC,1),size(PcritC,2),length(Rmf));
            d13_sp=nan(size(PcritC,1),size(PcritC,2),length(Rmf));
            for tt=1:length(Rmf)
            Rint(:,:,tt)=(Rmf(tt).*PcritC./pCO2ref+Rew)./(PcritC./pCO2ref+1);%length(Ac)*length(Eo)*length(Rmf)
            d13_int(:,:,tt)=(Rint(:,:,tt)./stn-1)*1e3;
            [~,~,~,~,Rarag,~]=C13speciation(Rint(:,:,tt),pHf,K1f,K2f,Env.Tfrac,spbC);
            d13_sp(:,:,tt)=(Rarag./stn-1)*1e3; 
            end
            
            %The isotopic values are now 3d if a vector of organic matter
            %endmembers is used
            D.Phic=Phic;D.PcritC=PcritC;D.VhPhic=VhPhic;D.PcritO2=PcritO2active;
            D.pco2w=pCO2ref;D.d13_dic=def;D.d13_int=d13_int;D.d13_sp=d13_sp;
            D.d13_aq=(Rew./stn-1)*1e3;D.d13_om=dmf-e_om;D.d13_met=dmf; 
            
            Qs(i,j,k)=qs;
                end %isotope model
                        
                % Average/variance across viable ecotypes
                AEloc=Trt.AEdist; %2d pdf of Ac and Eo
                AEloc(D.Phic<1 | D.Phic>Trt.phimax)=nan; 
                if length(Rmf)>1
                    WtTr=nan(size(AEloc,1),size(AEloc,2),length(Rmf));
                    for tt=1:length(Rmf)
                    WtTr(:,:,tt)=AEloc.*pD(tt); %trait weighting
                    end
                else
                    WtTr=AEloc;
                end
                
                d13o_wt=D.d13_sp.*WtTr; % otolith
                d13i_wt=D.d13_int.*WtTr; % internal ("blood")
                d13om_wt=D.d13_om.*pD; % organic matter without additional metabolic enrichment
                d13met_wt=D.d13_met.*pD; % metabolic endmember (if enrichment applied in settings)
                PcritO2_wt=D.PcritO2.*AEloc;
                PcritC_wt=D.PcritC.*AEloc;
                VhPhic_wt=D.VhPhic.*AEloc;
                Eo_wt=Trt.gE.*AEloc;
                
                %Note: this is average in each cell, not ocean average
                %(which requires volume weighting to avoid overrepresenting
                %surface and polar regions). The empircal distribution of d13oto
                %is saved as a distribution object for each hydrographic cell 

              % Empirical distribution of d13oto 
              clearvars biw Bw Cw kw ix iwi binW;
                [~,~,biw]=histcounts(D.d13_sp(:),bed); %take d13oto from any given cell
                Bw=nan(length(biw),1);
                %Redefine isotopic value to center of bin
                for iwi=1:length(biw)
                    if biw(iwi)~=0
                        Bw(iwi)=bcent(biw(iwi));
                    end
                end                
                %Generate weights in each bin
                Cw=[Bw WtTr(:)];
                kw=0; %initialize count through bins
                for iwi=bcent
                    ix=[];ix=Cw(:,1)==iwi;
                    kw=kw+1;
                    binW(kw)=sum(Cw(ix,2));   
                end
                wtPDFall(i,j,k,:)=binW;clearvars binW; %fractional weights for distb. in each location
          
              % Normal distribution parameters
                d13oto_3d_AEave(i,j,k)=nansum(d13o_wt(:))/nansum(WtTr(:));
                %d13oto_3d_AEstd(i,j,k)=std(d13o_wt(:)./WtTr(:),WtTr(:),'omitnan');
                d13int_3d_AEave(i,j,k)=nansum(d13i_wt(:))/nansum(WtTr(:));
                %d13int_3d_AEstd(i,j,k)=std(d13i_wt(:)./WtTr(:),WtTr(:),'omitnan');
                PcritO2_3d(i,j,k)=nansum(PcritO2_wt(:))/nansum(AEloc(:));
                PcritC_3d(i,j,k)=nansum(PcritC_wt(:))/nansum(AEloc(:));
                VhPhic_3d(i,j,k)=nansum(VhPhic_wt(:))/nansum(AEloc(:));
                Eo_3d(i,j,k)=nansum(Eo_wt(:))/nansum(AEloc(:));
                d13co2_3d(i,j,k)=D.d13_aq;
                d13dic_3d(i,j,k)=D.d13_dic;
                d13om_3d_DODave(i,j,k)=nansum(d13om_wt(:))/nansum(pD(:));
                %d13om_3d_DODstd(i,j,k)=std(d13om_wt(:)./pD(:),pD(:),'omitnan');
                d13met_3d_DODave(i,j,k)=nansum(d13met_wt(:))/nansum(pD(:));
                %d13met_3d_DODstd(i,j,k)=std(d13met_wt(:)./pD(:),pD(:),'omitnan');
                dic_3d(i,j,k)=Env.D;
                pco2_3d(i,j,k)=Env.C;   
                
               %Based on Anderson-Darling tests, d13oto is not normally
               %distributed in any cell. Thus, standard deviation not
               %meaningful and not retained (use empirical distribution).
            end
        end
    end
    toc;
end

%fix om and met so zeros replaced by nan
id=isnan(d13int_3d_AEave);
d13om_3d_DODave(id)=nan;
d13om_3d_DODstd(id)=nan;
d13met_3d_DODave(id)=nan;
d13met_3d_DODstd(id)=nan;
end

% Save outputs
if swSave==1
save d3D.mat d13oto_3d_AEave d13int_3d_AEave ...
    PcritO2_3d PcritC_3d VhPhic_3d d13co2_3d d13dic_3d d13om_3d_DODave...
    d13met_3d_DODave ...
    dic_3d pco2_3d WOAsub Qs Trt Eo_3d wtPDFall -v7.3
%d13oto_3d_AEstd d13om_3d_DODstd d13int_3d_AEstd d13met_3d_DODstd 

else
    load d3D.mat
end


