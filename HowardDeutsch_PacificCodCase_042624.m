%Author: Evan Howard (ehoward2@uw.edu)
%Original version: December 2020
%Current version: February 2023

%Description: This script compares temperature-dependent hypoxia based 
% models to observed Pacific cod otolith data, as a case study for
% time-resolved changes over lifespans.

%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 

%% Prepare workspace
clear all;close all;clc;
addpath('.../Helper_Scripts_and_files/'); %Location of helper scripts

%parameters to change
    %Bering Sea DIC & diet references: See Table S2
    %Juveniles may be like Atlantic cod, where they actually have the
    %heaviest dietary endmember because they are eating zooplankton that
    %are higher trophic level than the other diets (could be ~-18 per mil).
    %We don't make the trophic history that complicated for this toy model.
    ddic_larva=1.3;dom_larva=-22.5; % a little higher than lowest possible diet)
    ddic_adult=0.5;dom_adult=-17.5; % a little higher than lowest possible diet)

    %Assume a trophic enrichment of DIC produced from metabolism versus the
    %measured diet, tissue, or otolith organic matter?
    dOM=1; %change this as desired, but using same approximation as Atlantic cod
    dom_larva=dom_larva+dOM;dom_adult=dom_adult+dOM;
    
    %Equation constants
    kb=8.617e-5; Rgas=8.314e-3;T0=273.15;Tref=15;
    stn=0.011180; %VPDB 13C/12C

    %Estimated traits and endmembers
    Vh=0.05; %@whatever *environmental* reference temperature is used below
    %(this is not the same as Vh@15C)...this is a guess for 4C but could be
    %tuned
    Eo=0.35; %Gadus ogac (Greenland cod, may be same as Pacific cod 
    % according to some investigators) estimated as 0.25 @15C
    
    %%%
    dEdT=0.01; %temperature dependence of Eo (Atlantic cod -like value)
    ET=4; %temperature at which standard Eo is observed
    Eo=Eo+dEdT*(Tref-ET); % Eo value at Tref
    %%%
    
    Eovec=[0.25 0.5 0.7]; %arbitrary Eo at Tref
    Pc_med=2.5; %estimated PhiCrit for adult (something like interspecies median)
    PcL=1;PcA=Pc_med;
    pCO2_meas=400e-6;pO2ref=0.209;
    %Salinity in fish blood versus Bering Sea (as isionic correction)
    soc=32;sf=9; %summer salinity at 80 m depth--don't know seasonal cycle here
    pHb=7.85-0.15; %same approximation as in other analyses in this work

%Equation constants
kb=8.617e-5; Rgas=8.314e-3;T0=273.15;Tref=15;
stn=0.011180; %VPDB 13C/12C

%% Idealized time series
%Generate model ontogeny based on temperature and diet
%Ocean temperature
    %Adult tagging suggests that most of the sampled cod in Eastern Bering 
    %Sea and Gulf of Alaska were at 3-4C in March/April (low) to 6-7C in Sept/October (high), 
    %following a basically seasonal temperature pattern. This is not in 
    %agremeent with trawl catch temperatures, so the bottom temperatures 
    %from the trawl may not represent typical cod habitat. However, it 
    %could also be that adults regulate their environment (seasonal 
    %dependent DVM is observed, Nichol et al. 2013:10.1111/jfb.12160) so 
    %that larger surface temperature (0-9C with same seasonality) 
    %oscillations effect younger cod.

%Constant seasonal temperature cycle, adult
    T_octA=7;T_marA=3;%using maximum temperature cycle observed in adults
    funct_seasTA=@(x)mean([T_octA,T_marA])+(T_octA-T_marA)./2.*sin((2.*pi().*x+2.*pi.*9./12)); 
    %timeseries starts in March, phase shifted +3 months from origin
    %Constant seasonal temperature cycle, surface
    T_octS=9;T_marS=0;%hatch roughly in spring near low temperature
    funct_seasTS=@(x)mean([T_octS,T_marS])+(T_octS-T_marS)./2.*sin((2.*pi().*x+2.*pi.*9./12)); 
    %timeseries starts in March, phase shifted +3 months from origin
%Shifting seasonal cycle associated with lifestage
    t=[0:(1/12):7]';
    %guess for fraction of juvenile temperature & diet
    surt=8;trant=6; %8 months in surface, 6 month transition to depth
    fwt1=[ones(1,(surt-1)),linspace(1,0,(trant+1)),zeros(1,length(t)-((surt-1)+(trant+1)))]';
    %Make this a smooth logisitic function to avoid aesthetically displeasing
    %discontinuities
    funct_fwt=@(p,x)1./((1+exp(-p(1).*(x-p(2)))).^p(2)); pin_fwt=[-11.7;0.804];
    opts=statset('nlinfit');
    [param_fwt,~,~,~,~,~]=nlinfit(t,fwt1,funct_fwt,pin_fwt); 
    fwt=funct_fwt(param_fwt,t); 


%Ocean pH
%pH variability in the ocean has a minor effect compared to the pH
%difference between marine waters and fish blood, which is itself of
%secondary importance compared to most other temperature-dependent
%parameters in the otolith isotopic model. However, I include it for
%completeness. 
    %There is a surface pH record from the NOAA M2 buoy, which has 
    %April-June for a single year and May-October a second year. pH is 
    %highest in mid-May and lowest in mid-October (though no later data). 
    %Alternatively, use model output for 75 m depth, over shelf in E. Bering 
    %Sea--from Darren Pilcher's ROMS (pers. comm.). This has similar 
    %seasonal timing, but lower overall pH and a smaller amplitude, as 
    %might be expected at depth. For now using values estimated ROMS 
    % climatological model output for adult depths.
    %first month is March
    pHs=[7.84;7.91;7.93;7.97;7.97;7.89;7.82;7.75;7.72;7.69;7.7;7.73]; %Jan-Dec
    pHs=[pHs(3:12);pHs(1:2)]; %now time 0(i=1) is March
    funct_seaspH=@(x)pHs(round((x-floor(x)).*12+1));
    
%Ocean pCO2
%little carbonate system data available at depth, so using surface pCO2
%variability as stand (doi:10.1029/2005JC003074) but the same expectations 
%as with pH (above) probably apply. 
dpCO2=100e-6;
pCO2_octS=pCO2_meas-dpCO2;pCO2_marS=pCO2_meas+dpCO2;%hatch roughly in spring near low temperature
funct_seaspCO2S=@(x)mean([pCO2_octS,pCO2_marS])+(pCO2_octS-pCO2_marS)./2.*sin((2.*pi().*x+2.*pi.*9./12));
    
%Approximate fish mass with age
    %Von Bartalanffy equations and parameters for this Pacific Cod dataset 
    %(average across many fish, not specific to plotted individuals)
    load LMAparams;
    funct_W=@(x,p)p(1).*((1-exp(-p(2).*(x-p(3)))).^p(4));
    params=[LMAparams{4,1},LMAparams{2,1}(2),LMAparams{2,1}(3),LMAparams{1,1}(2)];
    w_t=funct_W(t,params); %This is only valid down to ~85 g mass, lower masses were not sampled
    iwt=w_t<85;
    w_t(w_t<85)=nan;w_t(1)=8e-5; %approximate mass at hatching: Laurel et al. 2008 (10.1093/plankt/fbn057)
    %use empirical, arbitrary fit to avoid negative numbers and smoothly
    %transition at lower mass
    funct_Ws=@(x,p)p(1).*(exp(-p(2)*x))+p(3);params_s=[-6.948,1.387,2.851]; %fit to log10 mass
    w10t=funct_Ws(t(iwt),params_s);w_t(iwt)=10.^w10t; %convert from log10 mass
    brange=w_t;
    %%%NOTE: ignoring biomass scaling of net O2 supply and demand (MI
    %Epsilon) because estimated changes over mass range are derived from
    %inter-species relationships that are not consistent with available
    %within-species otolith data. Available inter-species estimates for
    %various cod species are not far from Epsilon=0, suggesting little to
    %know size dependence beyond a couple hundred grams of mass. Since
    %there are no validated massess to pair with isotope composition
    %below 85 g, this cannot be calibrated in the current dataset and
    %assuming Epsilon=0 is the most conservative analytical choice.

%% Representative model
% First principles model of expected otolith evolution with examples
% showing that ontogenic diet change and seasonal temperature cycle are
% dominant factors explaining isotope evolution.

spbC='CO2aq'; %blood carbon species for mass balance 
%Assign reference conditions that are consistent with the approximate cold,
%high-biomass endmember 13C values in the Pacific cod otolith dataset
T_ref=4;B_ref=2.5e3;del_ref=-0.5; %vs T_ref=Tref=15C
    %%%Note: this is not the same as the default Tref used to calibrate to
    %%%measured hypoxia vulnerabilities at 15C in the MI database.
    %t(7)=4C, w_t(60)~2.5kg, d13Co~-0.5 is approximate asympototic d13Coto value   
brangeN=brange./B_ref;

RoP=cell(5,1);
for j=1:5
    if j==1 %only variation is pCO2 and pH
        %and allometric exponent (using adult values)
        clearvars trange Re Rm epsj Pc;
%         trange=funct_seasTA(t);
        trange=mean([T_octA,T_marA]).*ones(size(t));
        Re=(ddic_adult./1000+1).*stn;
        Rm=(dom_adult./1000+1).*stn;
        epsj=0;
        Pc=PcA;
    elseif j==2 %constant or variable allometric exponent over lifespan
        clearvars trange Re Rm epsj Pc;
        %trange=funct_seasTA(t);
        trange=mean([T_octA,T_marA]).*ones(size(t));
        Re=(ddic_adult./1000+1).*stn;
        Rm=(dom_adult./1000+1).*stn;
        Pc=PcA;
        epsj=-0.1; %nonzero epsilon
%     %Could calculate variable epsilon based on inter-species trend
%     load Allometry.mat; %Bg: centered log10 mass in g, delta, sigma, epsilon=delta-sigma (opposite convention used above and below)
%     B10c=abs(diff(Bg)./2)+Bg(2:end);
%     B10c=B10c(45:end)';epsilon=epsilon(45:end)';
%     funct_eps=@(p,x)p(1).*(1-1./(1+p(2).*exp(-x)).^p(3));pin_eps=[0.6126;2.0734;0.5701]; %fits log10w
%     opts=statset('nlinfit');opts.RobustWgtFun='logistic'; %Weight function used to minimize residual variances
%     [param_eps,~,~,~,~,~]=nlinfit(B10c,epsilon,funct_eps,pin_eps,opts); 
%     clearvars epsilon delta Bc B10c opts; 
%         epsj=-1.*funct_eps(param_eps,log10(brange));       
    elseif j==3 %variable T seasonality over lifespan   
        clearvars trange Re Rm epsj Pc;
        trange=fwt.*funct_seasTS(t)+(1-fwt).*funct_seasTA(t);
        Re=(ddic_adult./1000+1).*stn; %guess Re
        Rm=(dom_adult./1000+1).*stn; %guess Rm
        epsj=0;
        Pc=PcA;
    elseif j==4 %variable isotopic endmembers over lifespan
        clearvars trange Re Rm epsj Pc;
        %trange=funct_seasTA(t);
        trange=mean([T_octA,T_marA]).*ones(size(t));
        Re=fwt.*((ddic_larva./1000+1).*stn)+(1-fwt).*((ddic_adult./1000+1).*stn); %guess Re
        Rm=fwt.*((dom_larva./1000+1).*stn)+(1-fwt).*((dom_adult./1000+1).*stn); %guess Rm
        epsj=0;
        Pc=PcA;
    elseif j==5 %variable T, endmembers, and epsilon (full model)
        clearvars trange Re Rm epsj;
        trange=fwt.*funct_seasTS(t)+(1-fwt).*funct_seasTA(t);
        Re=fwt.*((ddic_larva./1000+1).*stn)+(1-fwt).*((ddic_adult./1000+1).*stn); %guess Re
        Rm=fwt.*((dom_larva./1000+1).*stn)+(1-fwt).*((dom_adult./1000+1).*stn); %guess Rm
        epsj=-0.1;
        %epsj=-1.*funct_eps(param_eps,log10(brange)); %If variable Epsilon
        %over time is modeled instead
        Pc=PcA;
        %Pc=fwt.*PcL+(1-fwt).*PcA; %from 1 to 2.5, as example
        %Can change active to resting ratio (PhiCrit), e.g., could shift 
        % with age or temperature. However, based on sensitivity analyses
        %for this case study, either this isn't sufficient to drive key 
        % features of time series, and there aren't independent constraints. 
        % Something to investigate in the future.      
    end

    TT=(1/kb).*(1./(trange+T0)-1./(T_ref+T0)); %temperature scaling 
    % relative to Vh at environmental reference temperature
    
%Aqueous chemistry
%What is the new pH because of temperature change only?
pHe=pH_T(trange,soc,funct_seaspH(t),trange);pHe=pHe';
pHf=pH_T(trange,sf,pHb,trange);pHf=pHf';
if (j~=1) & (j~=5)
    pHe=repmat(mean(pHe),length(trange),1);pHf=repmat(mean(pHf),length(trange),1);
end
%carbon dioxide solubility
[K0e K1e K2e]=carbeq(trange,soc);K0e=K0e';K1e=K1e';K2e=K2e'; %In seawater
[K0f K1f K2f]=carbeq(trange,sf);K0f=K0f';K1f=K1f';K2f=K2f'; %In fish bodily fluids

%carbonate system partioning
He=10.^(-1.*pHe);De=(He.^2)+He.*K1e'+K1e'.*K2e';
X0e=(He.^2)./De;X1e=He.*K1e'./De;%X2=1-X0-X1;
X2e=K1e'.*K2e'./De;
Hf=10.^(-1.*pHf);Df=(Hf.^2)+Hf.*K1f'+K1f'.*K2f';
X0f=(Hf.^2)./Df;X1f=Hf.*K1f'./Df;%X2=1-X0-X1;
X2f=K1f'.*K2f'./Df;
%oxygen solubility
Khe=O2sol(soc,trange)./(pO2ref).*1e-6;Khe=Khe'; %mol kg-1 atm-1 of O2
Khf=O2sol(sf,trange)./(pO2ref).*1e-6;Khf=Khf'; 
%Diffusivities in water (all in 10-4 m2 s-1)
[Dco2e,Dhco3e,Dco3e]=co2_diffusion(soc,trange);Dco2e=Dco2e';Dhco3e=Dhco3e';Dco3e=Dco3e';
[Dco2f,Dhco3f,Dco3f]=co2_diffusion(sf,trange);Dco2f=Dco2f';Dhco3f=Dhco3f';Dco3f=Dco3f';
[~,~,~,~,Do2e,~,~]=gas_diffusion(soc,trange);Do2e=Do2e';
[~,~,~,~,Do2f,~,~]=gas_diffusion(sf,trange); Do2f=Do2f';

%isotopic endmembers
tmpRe=nan(length(trange),1);
tmpRm=nan(length(trange),1);
for i=1:length(trange)
    if length(Re)==1
    [~,tReG,~,~,~,~]=C13speciation(Re,pHe(i),K1e(i),K2e(i),trange(i),'DIC'); %CO2aq,e   
    tmpRe(i)=tReG;tmpRm(i)=Rm;clearvars tReG;
    else
    [~,tReG,~,~,~,~]=C13speciation(Re(i),pHe(i),K1e(i),K2e(i),trange(i),'DIC'); %CO2aq,e   
    tmpRe(i)=tReG;tmpRm(i)=Rm(i);clearvars tReG; 
    end
end

%stoichiometry
qd=1;
qs=K0e.*Dco2e./(Khe.*Do2e);%Vh new vs Vh old (blood pCO2 effect on estimate)
fs=1;fd=0.5; %factors to convert original Vh to blood dependence
qs=qs.*fs;qd=qd.*fd;
qsT=qs';   
%Calculate otolith composition
clearvars RoCO2;
tt=TT;p_m=tmpRm(:);p_e=tmpRe(:);
Eotmp=Eo+dEdT*(Tref-15); % Eo value at new Tref (4C) relative to prior at 15
em=Eotmp+dEdT.*(trange-Tref);
pCO2tmp=funct_seaspCO2S(t);
if (j~=1) & (j~=5)
    pCO2tmp=mean(pCO2tmp);
end

RoCO2=((p_m.*Vh.*Pc./pCO2tmp.*qd./qsT.*exp(-em.*tt)).*brangeN.^(-epsj)+p_e)./...
        ((Vh.*Pc./pCO2tmp.*qd./qsT.*exp(-em.*tt)).*brangeN.^(-epsj)+1);  
for i=1:length(trange)
[~,~,~,~,R13Arag,~]=C13speciation(RoCO2(i),pHf(i),K1f(i),K2f(i),trange(i),spbC);
RoP{j}(i,1)=(R13Arag./stn-1)*1e3;
end
clearvars tt p_m p_e allRm allRe TTb i; 
end
%% Plot model and data
%Plot model
%Don't know how much 'preformed' (adult composition) carbon is in the 
%larval fish and averaged with the larval fishes' contributions during 
%milling; the samples had a width increment about 10 times larger in the 
%otolith core than when trenching to the exterior edge, so first increment 
%is definitely averaging over a lot of material that may have composition trend.
%To focus on where model could actually have data constraints, just trim
%first six months off model timeseries.
it=7:length(t);
iw=find(brange>=85);%the lowest mass measured was ~85 g in this dataset. 

%figure graphical properties
figure('Renderer', 'painters', 'Position', [100 100 250 750]);
subplot(313); %Model outputs
%Only carbonate system seasonality
plot(t(it),RoP{1}(it),'-','Color','b','LineWidth',1.5);hold on;
%Changing epsilon over lifespan
plot(t(it),RoP{2}(it),'-','Color',[0.6 0 1],'LineWidth',1.5);
%Changing T seasonality over lifespan
plot(t(it),RoP{3}(it),'-','Color','r','LineWidth',1.5);
 %Changing diet over lifespan
plot(t(it),RoP{4}(it),'-','Color',[1 0.5 0],'LineWidth',1.5);
%Full lifespan transition
plot(t(it),RoP{5}(it),'-','Color','k','LineWidth',2);
xlabel('Model age (y)');ylabel('\delta^{13}C_{o}');
ylim([-6 1]);xlim([0.5 3.5]);

subplot(312); %Model forcings
plot(t(it),trange(it),'-','Color','r','LineWidth',1.5);hold on;%temperature
plot(t(it),log10(brange(it)),'Color',[0.6 0 1],'LineWidth',1.5);%biomass
ylabel('\color[rgb]{0.6 0 1}log_{10}Mass (g), \color{red}T (^{o}C)');xlabel('Model age (y)');
yyaxis right;
plot(t(it),(tmpRm(it)./stn-1)*1e3,'-','Color',[1 0.5 0],'LineWidth',1.5);%diet
ylabel('\delta^{13}C_{met}');
xlim([0.5 3.5]);
ax=gca;ax.YAxis(2).Color=[1 0.5 0];ax.YAxis(1).Color='r';

subplot(311); %example data
%Pacific cod data from Kastelle et al.
load HaulInfo;%Cruise,Haul,Specimen,Sex(unk),Length(mm),Weight(g),T_bottom,T_surface
load IsotopeInfo; %Cruise,Haul,Specimen,winter #, width of band (mm), 
                  %running sum of width from initial track (mm), d13C, d18O
%correct various transcription errors in data table
IsotopeInfo(33,4)=5;IsotopeInfo(99,4)=3;IsotopeInfo(295,4)=2;
IsotopeInfo(360,4)=4;IsotopeInfo(428,4)=4;IsotopeInfo(509,4)=5;
IsotopeInfo(1055,4)=4;  %fix years for entries with nonstandard notation
IsotopeInfo(204,5)=0.33; %fix incorrect width
IsotopeInfo([511,1141],:)=[]; %remove missing rows
%find unique samples
IDh=string(HaulInfo(:,1:3));IDh=append(IDh(:,1),IDh(:,2),IDh(:,3));
ID=string(IsotopeInfo(:,1:3));ID=append(ID(:,1),ID(:,2),ID(:,3));
IDu=unique(ID,'stable');
%subset to just a few examples to illustrate isotopic evolution in otolith
IDu=IDu([7,16,22,36]);
ij=zeros(length(ID),1);ih=zeros(length(IDh),1);
for i=1:length(IDu)
ii=find(ID==IDu(i));ij(ii)=1;clearvars ii;
ii=find(IDh==IDu(i));ih(ii)=1;
end
%remove records and increments with missing data
ij=logical(ij);ih=logical(ih);IsotopeInfo(~ij,:)=[];ID(~ij)=[];HaulInfo(~ih,:)=[];IDh(~ih)=[];
clearvars ii ij ih i;
ID(isnan(IsotopeInfo(:,7)))=[];IsotopeInfo(isnan(IsotopeInfo(:,7)),:)=[];
%identify first winter peak based on oxygen isotope record, as identified
%in Kastelle et al.
ip=cell(4,1);
ip{1}(1,1)=10;ip{2}(1,1)=10;ip{3}(1,1)=10;ip{4}(1,1)=7; %ip{3} could be 3
%subtract relative otolith distance at winter 1 to normalize lengths
apxWidth=IsotopeInfo(:,6);apxWidth(:,2)=nan;
for i=1:length(IDu)
    k=find(ID==IDu(i));tmp=IsotopeInfo(k,6);
    %center on first winter, don't try to do other stretching as don't
    %know individual growth curves
    apxWidth(k,2)=(apxWidth(k,1)-tmp(ip{i}(1)));
end
clearvars tmp k i;
%make the plot
Clr=[0.7 0.1 0.7;1 0 0;1 0.7 0;0 0 1];
for i=1:length(IDu)
    k=find(ID==IDu(i));
    tx=apxWidth(k,2);
    ty=IsotopeInfo(k,7);
    clr=Clr(i,:);
    plot(tx,ty,'.-','Color',clr,'MarkerEdgeColor',clr,'LineWidth',1.5,'MarkerSize',12);hold on;
end
clearvars i k tx ty tx2 ty2 ix2
xlabel('\Delta width (mm)');ylabel('\delta^{13}C_{o}');    
xlim([-0.75 1.5]);ylim([-6 1]);
