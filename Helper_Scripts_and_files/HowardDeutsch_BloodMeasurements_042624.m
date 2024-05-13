%Author: Evan Howard (ehoward2@uw.edu)
%Original version: February 2022
%Current version: August 2022

%Description: This script corrects reported blood and extracellular pCO2 
%(Pint) for differences in conditions and carbonate system assumptions across
%published experiments, and plots the results adjusted to 15C as groups
%related to resting and excercised status and the point of circulation
%sampled (e.g., blood from pre- or post-gills). Only pCO2 from teleost fish 
%or closely related Actinopterygii are included, though blood pCO2 from 
%elasmobranchs and agnathates seems broadly similar.

%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 

%% Prepare workspace
clear all;
swW=0; %if 1, write in literature reported values
%% Write data
if swW==1
%General Note: In many of these references pCO2 was calculated from unreported 
%measurements of total inorganic carbon using the solubility and speciation
%relationships described in Boutilier et al. 1984. However, these
%relationships are not significantly different than those for seawater when
%corrected for temperature and ionic strength differences, and the latter
%are more accurately known.

%If reference does not provide pCO2 of water in experiment, but 
%specifically notes aeration, efforts to maintain atmospheric 
%equilibrium, etc, set reported pCO2 to 400 uatm--slightly elevated from
%atmospheric for most of these periods--because measured pCO2 in
%respirometry experiments is ALWAYS supersaturated (no one controls for
%constant carbonate chemistry in water). If any indication of less than
%optimal aeration approach, set to 800 uatm, typical of available water
%measurements for respirometry. If used tap water, recirculated water, or
%other source with reported depletion in pH indicating high in line
%respiration prior to experiment, set to 2000 uatm, typical of municipal
%tap water sources (see inline references below).

%The difference between 300-500 uatm is not really important to the 
%calculations here, the difference between near atmospheric and 
%substantially elevated pCO2 (typical experimental conditions) strongly 
%influences estimated Pmet/Pw (potential biases of factor of 2-5), and 
%weakly influences estimated active/resting ratio (order of 20%).

for i=1:100
    if i==1
Reference="Lee 2003";
Species="Seriola quinqueradiata";
Environment="Marine";
Fluid="Low-Dorsal aorta";
Activity="Low swim speed";
Measurement="Calculated";
%Note: Using initial point only, pre-acidification experiment
T_rep=20;
P_rep=3/760*1e6; %mmHg to uatm, calculated
P_w_rep=800; %"well aerated seawater" at pH=8.25, but actual experiment in
%closed/not aerated respirometry chambers, so 800uatm more plausible than
%400 uatm.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end
    
    if i==2
Reference="Larsen 1997";
Species="Gadus morhua";
Environment="Marine";
Fluid="High-Ventral aorta";
Measurement="Measured";
Activity="Rest";
%Note: Using initial point only, pre-acidification experiment
T_rep=12;
P_rep=3.5/760*1e6; %mmHg to uatm, measured
P_w_rep=400; %"well-aerated seawater"
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end
    
    if i==3
Reference="McKenzie 2002";
Species="Anguilla anguilla";
Environment="Fresh";
Fluid="Low-Dorsal aorta";
Measurement="Calculated";
Activity="Rest";
%Note: Using initial point only, pre-acidification experiment
T_rep=23; %But TIC only measured at 37C; assuming here T-corrected to 23, but not clear
P_rep=3.0/760*1e6; %mmHg to uatm
P_w_rep=0.6/760*1e6; %measured
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end

    if i==4
Reference="Michaelidis 2007";
Species="Sparus aurata";
Environment="Marine";
Fluid="Mod-Caudal vein";
Measurement="Calculated";
Activity="Starved";
%Note: Using initial point only, pre-acidification experiment
T_rep=18; %Fish at 15C, but pH measured at 18C and TIC at 37C. Paper not
%sufficiently clear to know what temperatures different elements were
%corrected to; I assume here that TIC was corrected to 18C, becuase no
%indication of where or how pH would have been temperature corrected.
P_rep=2.56/760*1e6; %mmHg to uatm
P_w_rep=0.82/760*1e6; %measured
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end

    if i==5
Reference="Thomas 1983";
Species="Oncoryhnchus mykiss";
Environment="Fresh";
Fluid="Low-Subclavian artery";
Measurement="Measured";
Activity="Low swim speed";
%Note: Using initial point only, pre-acidification experiment
T_rep=15; %Fish at 15C, 
P_rep=2.3/760*1e6; %mmHg to uatm, measured
P_w_rep=0.25/760*1e6; %measured
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end

    if i==6
Reference="Thomas 1987";
Species="Oncoryhnchus mykiss";
Environment="Fresh";
Fluid="Low-Subclavian artery";
Measurement="Measured";
Activity="High swim speed";
%Lots going on here, so just selecting the 50cm (Ucrit) fit on Fig. 2
%Note: Using initial point only, pre-acidification experiment
T_rep=12;
P_rep=(2+1)/760*1e6; %mmHg to uatm, measured
P_w_rep=0.67/760*1e6; %measured
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end
    
    if i==7
Reference="Bernier 2004";
Species="Oncoryhnchus tshawytscha";
Environment="Marine";
Fluid="Low-Dorsal aorta";
Measurement="Calculated";
Activity="High swim speed";
%Using high swim speed values, e.g., Ucrit
T_rep=9;
P_rep=(4.3)/760*1e6; %mmHg to uatm, calculated
P_w_rep=800; %kept in low-flow chamber respirometer for day prior to experiments
%Thus likely to have had high pCO2 of environmental water. Similar
%experiments in this and other respirometry datasets have average pCO2
%order of 800-1000 uatm, or two to three times atmospheric.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end    
    
    if i==8
Reference="Bernier 2004";
Species="Oncoryhnchus tshawytscha";
Environment="Marine";
Fluid="Low-Dorsal aorta";
Measurement="Calculated";
Activity="Low swim speed";
%Using low swim speed values, e.g., 0.4 body lengths s-1
T_rep=9;
P_rep=(3)/760*1e6; %mmHg to uatm, calculated
P_w_rep=800; %kept in low-flow chamber respirometer for day prior to experiments
%Thus likely to have had high pCO2 of environmental water. Similar
%experiments in this and other respirometry datasets have average pCO2
%order of 800-1000 uatm, or two to three times atmospheric.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end  
    
    if i==9
Reference="Wood and Eom 2019";
Species="Oncoryhnchus mykiss";
Environment="Fresh";
Fluid="Low-Dorsal aorta";
Measurement="Measured";
Activity="Ventilated/Sedated";
%Using starved fish dataset since ventilated, not swimming anyway
T_rep=12;
P_rep=(3.6)/760*1e6; %mmHg to uatm, measured
P_w_rep=2000; %Vancouver BC tap water of pH=7, whereas reported Vancouver 
%tap water is buffered to pH=7.5 by the city; thus likely significant 
%respiration or other causes of very high pCO2 in water. Non-tropical river/tap water
%with similar decreases in pH corresponds to pCO2 on order of >2000uatm
%in other settings (e.g., https://bg.copernicus.org/articles/12/67/2015/).
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end  
    
    if i==10
Reference="Wood and Eom 2019";
Species="Carassius auratus";
Environment="Fresh";
Fluid="Low-Dorsal aorta";
Measurement="Measured";
Activity="Ventilated/Sedated";
%Using starved fish dataset since ventilated, not swimming anyway
T_rep=12;
P_rep=(5)/760*1e6; %mmHg to uatm, measured
P_w_rep=2000; %Vancouver BC tap water of pH=7, whereas reported Vancouver 
%tap water is buffered to pH=7.5 by the city; thus likely significant 
%respiration or other causes of very high pCO2 in water. Non-tropical river/tap water
%with similar decreases in pH corresponds to pCO2 on order of >2000uatm
%in other settings (e.g., https://bg.copernicus.org/articles/12/67/2015/, 
%https://doi.org/10.3390/rs13234916)
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end    
    
    if i==11
Reference="Gam 2020";
Species="Chitala ornata";
Environment="Fresh";
Fluid="Low-Dorsal aorta";
Measurement="Measured";
Activity="Rest";
%Using only lowest T, since appreciable air-breathing (gulps not gills) at higher T
T_rep=25;
P_rep=(5)/760*1e6; %mmHg to uatm, measured, but some imprecision in how reported
P_w_rep=900; %"constant aeration, pCO2<0.7mmHg, pH 7.75-7.85". So<920 uatm.
%Rounding to 900 since 920 is max but water pH indicates can't be much
%lower.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end     
    
    if i==12
Reference="Ern 2016";
Species="Sciaenops ocellatus";
Environment="Marine";
Fluid="Mod-Caudal artery";
Measurement="Calculated";
Activity="Rest";
%Using average of two 380 uatm controls in Table 1
T_rep=22;
P_rep=(3.08+3.18)/2/760*1e6; %uatm, calculated
P_w_rep=380; %reported based on design, but not explicitly measured
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end       
    
    if i==13
Reference="Kwan 2022";
Species="Sebastes diploproa";
Environment="Marine";
Fluid="Mod-Caudal artery";
Measurement="Calculated";
Activity="Ventilated/Sedated";
%Using lower pCO2 treatment
T_rep=18;
P_rep=1603; %uatm, calculated
HCO3_rep=2.37; %mmol L-1, TIC measured, HCO3 calculated and reported
pH_rep=7.75; %measured, NBS scale
P_w_rep=572; %measured
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end         
    
    if i==14
Reference="Montgomery 2022";
Species="Dicentrarchus labrax";
Environment="Marine";
Fluid="Mod-Caudal artery";
Measurement="Calculated";
Activity="Rest";
%Using control conditions
T_rep=13.94;
P_rep=0.2/101.325*1e6; %uatm, calculated
P_w_rep=0.059/101.325*1e6; %kPa to uatm, measured
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end     
    
    if i==15
Reference="Hayashi 2004";
Species="Paralichthys olivaceus";
Environment="Marine";
Fluid="Mod-Caudal artery";
Measurement="Calculated";
Activity="Starved";
%Using control conditions for flounder only, data for other two species 
%not fully reported (and anyway,  one is elasmobranch and other may be 
%duplicate data with reference #1).
T_rep=20;
P_rep=0.25/101.325*1e6; %kPa to uatm, calculated
P_w_rep=400; %Not reported, "well aerated seawater", "fresh seawater into experimental tank"
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end 
    
    if i==16
Reference="Toews 1983";
Species="Conger conger";
Environment="Marine";
Fluid="Low-Dorsal aorta";
Measurement="Calculated";
Activity="Starved";
T_rep=17;
P_rep=2/760*1e6; %mmHg to uatm, calculated
P_w_rep=0.6/760*1e6; %"0.5 to 0.7mmHg"
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end    
    
%NOTE: i=17:20 are low/high swim speed, venous/arterial blood from 
%same experiments.
    if i==17
Reference="Korsmeyer 1997a and 1997b";
Species="Thunnus albacares"; %CAUTION: warm-blooded tuna, may not be 
%be suitable for otolith model predictions/comparisons. Experimental
%conditions in 1997a, blood details in 1997b
Environment="Marine";
Fluid="Low-Dorsal aorta";
Measurement="Measured";
Activity="Slow swim";
T_rep=25; %fish T 24C, but all measurements at 25C
P_rep=0.81/101.325*1e6; %kPa to uatm, calculated
P_w_rep=800; %not reported, recirculated warm water so likely high pCO2
%estimating as ~800 uatm similar to other swim tunnel experiments
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end     
    
     if i==18
Reference="Korsmeyer 1997a and 1997b";
Species="Thunnus albacares"; %CAUTION: warm-blooded tuna, may not be 
%be suitable for otolith model predictions/comparisons. Experimental
%conditions in 1997a, blood details in 1997b
Environment="Marine";
Fluid="Low-Dorsal aorta";
Measurement="Measured";
Activity="Fast swim";
T_rep=25; %fish T 24C, but all measurements at 25C
P_rep=0.89/101.325*1e6; %kPa to uatm, calculated
P_w_rep=800; %not reported, recirculated warm water so likely high pCO2
%estimating as ~800 uatm similar to other swim tunnel experiments
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
     end    
    
     if i==19
Reference="Korsmeyer 1997a and 1997b";
Species="Thunnus albacares"; %CAUTION: warm-blooded tuna, may not be 
%be suitable for otolith model predictions/comparisons. Experimental
%conditions in 1997a, blood details in 1997b
Environment="Marine";
Fluid="High-Ventral aorta";
Measurement="Measured";
Activity="Slow swim";
T_rep=25; %fish T 24C, but all measurements at 25C
P_rep=0.91/101.325*1e6; %kPa to uatm, calculated
P_w_rep=800; %not reported, recirculated warm water so likely high pCO2
%estimating as ~800 uatm similar to other swim tunnel experiments
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
     end  
    
     if i==20
Reference="Korsmeyer 1997a and 1997b";
Species="Thunnus albacares"; %CAUTION: warm-blooded tuna, may not be 
%be suitable for otolith model predictions/comparisons. Experimental
%conditions in 1997a, blood details in 1997b
Environment="Marine";
Fluid="High-Ventral aorta";
Measurement="Measured";
Activity="Fast swim";
T_rep=25; %fish T 24C, but all measurements at 25C
P_rep=1.1/101.325*1e6; %kPa to uatm, calculated
P_w_rep=800; %not reported, recirculated warm water so likely high pCO2
%estimating as ~800 uatm similar to other swim tunnel experiments
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
     end   
    
    if i==21
Reference="van den Thillart 1983";
Species="Oncorhynchus kisutch";
Environment="Marine"; %Caution, some experiments with manipulated water 
%chemistry, not using those.
Fluid="Low-Dorsal aorta";
Measurement="Measured";
Activity="Slow swim"; %"Rest"
T_rep=13;
P_rep=2.5/760*1e6; %mmHg to uatm, measured
P_w_rep=800; %lots of weird buffered chemistry, details for the 
%normal seawater trials are insufficent to evaluate PCO2w. This is a
%naturally high pCO2 area (>600 uatm would be typical) and in respirometer,
%so assume 800 uatm as with other unknown respirometry trials.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end    
    
    if i==22
Reference="van den Thillart 1983";
Species="Oncorhynchus kisutch";
Environment="Marine"; %Caution, some experiments with manipulated water 
%chemistry, not using those.
Fluid="Low-Dorsal aorta";
Measurement="Measured";
Activity="Fast swim"; %"Burst swim", above Ucrit
T_rep=13;
P_rep=8.0/760*1e6; %mmHg to uatm, measured
P_w_rep=800; %lots of weird buffered chemistry, details for the 
%normal seawater trials are insufficent to evaluate PCO2w. This is a
%naturally high pCO2 area (>600 uatm would be typical) and in respirometer,
%so assume 800 uatm as with other unknown respirometry trials.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end       
    
    if i==23
Reference="Holeton 1983";
Species="Oncoryhnchus mykiss";
Environment="Fresh";
Fluid="Low-Dorsal aorta";
Measurement="Measured";
Activity="Slow swim"; %initial condition, compared to next which was fast swim
%stimulated by electric schocks (these poor fish)
T_rep=15;
P_rep=2/760*1e6; %mmHg to uatm, measured
P_w_rep=0.6/760*1e6; %"<0.6 mm Hg" in respirometer setup (higher than holding tanks)
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end      
    
    if i==24
Reference="Holeton 1983";
Species="Oncoryhnchus mykiss";
Environment="Fresh";
Fluid="Low-Dorsal aorta";
Measurement="Measured";
Activity="Fast swim"; %T0 after stimulated by electric schocks (these poor fish)
T_rep=15;
P_rep=3/760*1e6; %mmHg to uatm, measured
P_w_rep=0.6/760*1e6; %"<0.6 mm Hg" in respirometer setup (higher than holding tanks)
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end  
    
    if i==25
Reference="Milligan 1987";
Species="Oncoryhnchus mykiss";
Environment="Fresh";
Fluid="Low-Dorsal aorta";
Measurement="Calculated";
Activity="Rest"; %using initial point of excercise group, since 'control' group
%started in slightly different place rather than a carbonate system control
T_rep=15; %sample bath temperature, close to fish expt temperature
P_rep=2.7/760*1e6; %mmHg to uatm, calculated
P_w_rep=800; %Again tap water with low pH, but says "well-aerated", 
%so I'm assuming 800 uatm rather than near-atmosphere or the the much
%higher values for other tap waters; e.g. intermediate between initial
%tapwater and atmospheric equilibrium. This is probably reasonable because
%well aerated tap waters with respect to O2 still tend to be supersaturated
%with respect CO2.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end  
    
    if i==26
Reference="Milligan 1987";
Species="Oncoryhnchus mykiss";
Environment="Fresh";
Fluid="Low-Dorsal aorta";
Measurement="Calculated";
Activity="Fast swim"; %T0 after chased
T_rep=15; %sample bath temperature, close to fish expt temperature
P_rep=5.5/760*1e6; %mmHg to uatm, calculated
P_w_rep=800; %Again tap water with low pH, but says "well-aerated", 
%so I'm assuming 800 uatm rather than near-atmosphere or the the much
%higher values for other tap waters; e.g. intermediate between initial
%tapwater and atmospheric equilibrium. This is probably reasonable because
%well aerated tap waters with respect to O2 still tend to be supersaturated
%with respect CO2.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end     
    
    if i==27
Reference="Milligan 1987";
Species="Platichthys stellatus";
Environment="Marine";
Fluid="Mod-Caudal artery";
Measurement="Calculated";
Activity="Rest"; %Again using excercise group Tinitial because group differences
T_rep=10; %sample bath temperature, close to fish expt temperature of 9
P_rep=2.4/760*1e6; %mmHg to uatm, calculated
HCO3_rep=6.7; %mmol L-1, HCO3 calculated from TIC (TIC not reported)
pH_rep=7.79; %measured, NBS scale assumed but not stated
P_w_rep=800; %Not reported: Friday Harbor Labs seawater from late fall is order of 800 uatm; 
%doi: 10.1002/lno.10062. But could be even higher because in tank setup
%without flow, so tank gas exchange may not be fast compared to respiration
%in these experiments.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;

    end  
    
    if i==28
Reference="Milligan 1987";
Species="Platichthys stellatus";
Environment="Marine";
Fluid="Mod-Caudal artery";
Measurement="Calculated";
Activity="Active"; %Excercise T=0 min post excercise
T_rep=10; %sample bath temperature, close to fish expt temperature of 9
P_rep=5.5/760*1e6; %mmHg to uatm, calculated
HCO3_rep=7.2; %mmol L-1, HCO3 calculated from TIC (TIC not reported)
pH_rep=5.5; %measured, NBS scale assumed but not stated
P_w_rep=800; %Not reported: Friday Harbor Labs seawater from late fall is order of 800 uatm; 
%doi: 10.1002/lno.10062. But could be even higher because in tank setup
%without flow, so tank gas exchange may not be fast compared to respiration
%in these experiments.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end  
    
    if i==29 
Reference="Currie and Tufts 1993";
Species="Oncorhynchus mykiss";
Environment="Fresh";
Fluid="Low-Dorsal artery";
Measurement="Calculated";
Activity="Resting"; 
T_rep=10; %fish at 9-11
P_rep=0.31/101.325*1e6; %kPa to uatm, calculated
P_w_rep=800; %Not reported, "aerated fresh water". Given context likely tap water as well,
%so 800-2000 uatm most likely
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end
    
    if i==30 
Reference="Currie and Tufts 1993";
Species="Oncorhynchus mykiss";
Environment="Fresh";
Fluid="Low-Dorsal artery";
Measurement="Calculated";
Activity="Active"; 
T_rep=10; %fish at 9-11
P_rep=0.67/101.325*1e6; %kPa to uatm, calculated
P_w_rep=800; %Not reported, "aerated fresh water". Given context likely tap water as well,
%so 800-2000 uatm most likely
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end    
    
    if i==31 
Reference="Currie and Tufts 1993";
Species="Oncorhynchus mykiss";
Environment="Fresh";
Fluid="High-Ventral aorta";
Measurement="Calculated";
Activity="Resting"; 
T_rep=10; %fish at 9-11
P_rep=0.41/101.325*1e6; %kPa to uatm, calculated
P_w_rep=800; %Not reported, "aerated fresh water". Given context likely tap water as well,
%so 800-2000 uatm most likely
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end      
    
    if i==32 
Reference="Currie and Tufts 1993";
Species="Oncorhynchus mykiss";
Environment="Fresh";
Fluid="High-Ventral aorta";
Measurement="Calculated";
Activity="Active"; 
T_rep=10; %fish at 9-11
P_rep=1.15/101.325*1e6; %kPa to uatm, calculated
P_w_rep=800; %Not reported, "aerated fresh water". Given context likely tap water as well,
%so 800-2000 uatm most likely
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end   
    
    if i==33 
Reference="Egginton 1997"; %Using only the two fish with Hemoglobin, 
%not trying to generalize this to icefish with no blood cells yet.
%Care/water details from Egginton 1994.
Species="Notothenia coriiceps";
Environment="Marine";
Fluid="Low-Dorsal aorta";
Measurement="Calculated";
Activity="Resting"; 
T_rep=1; %unclear what fish were at, but chemistry at 1C
P_rep=0.35/101.325*1e6; %kPa to uatm, calculated
P_w_rep=800; %not reported. Buckets on deck at ~1C. Tank for excercise not described. 
%my assumption is that pCO2 was likely to be elevated compared to
%ambient Southern ocean water (~400 for that time of year, a little 
%supersaturated wrt atmosphere). Given lack of flow and reused water, 800
%uatm seems reasonable, probably not much higher but possibly lower.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end 
    
    if i==34 
Reference="Egginton 1997"; %Using only the two fish with Hemoglobin, 
%not trying to generalize this to icefish with no blood cells!
%Holding details from Egginton 1994.
Species="Notothenia coriiceps";
Environment="Marine";
Fluid="Low-Dorsal aorta";
Measurement="Calculated";
Activity="Active"; 
T_rep=1; %unclear what fish were at, but chemistry at 1C
P_rep=0.7/101.325*1e6; %kPa to uatm, calculated
P_w_rep=800; %not reported. Buckets on deck at ~1C. Tank for excercise not described. 
%my assumption is that pCO2 was likely to be elevated compared to
%ambient Southern ocean water (~400 for that time of year, a little 
%supersaturated wrt atmosphere). Given lack of flow and reused water, 800
%uatm seems reasonable, probably not much higher but possibly lower.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end  
    
    if i==35 
Reference="Egginton 1997"; %Using only the two fish with Hemoglobin, 
%icefish from same paper not used as unclear whether model applies (but would
%nonetheless have a similar activity ratio to these two rockfish).
%Holding details from Egginton 1994.
Species="Notothenia rossii";
Environment="Marine";
Fluid="Low-Dorsal aorta";
Measurement="Calculated";
Activity="Resting"; 
T_rep=1; %unclear what fish were at, but chemistry at 1C
P_rep=0.35/101.325*1e6; %kPa to uatm, calculated
P_w_rep=800; %not reported. Buckets on deck at ~1C. Tank for excercise not described. 
%my assumption is that pCO2 was likely to be elevated compared to
%ambient Southern ocean water (~400 for that time of year, a little 
%supersaturated wrt atmosphere). Given lack of flow and reused water, 800
%uatm seems reasonable, probably not much higher but possibly lower.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end   
    
    if i==36 
Reference="Egginton 1997"; %Using only the two fish with Hemoglobin, 
%not trying to generalize this to icefish with no blood cells!
%Holding details from Egginton 1994.
Species="Notothenia rossii";
Environment="Marine";
Fluid="Low-Dorsal aorta";
Measurement="Calculated";
Activity="Active"; 
T_rep=1; %unclear what fish were at, but chemistry at 1C
P_rep=0.55/101.325*1e6; %kPa to uatm, calculated
P_w_rep=800; %not reported. Buckets on deck at ~1C. Tank for excercise not described. 
%my assumption is that pCO2 was likely to be elevated compared to
%ambient Southern ocean water (~400 for that time of year, a little 
%supersaturated wrt atmosphere). Given lack of flow and reused water, 800
%uatm seems reasonable, probably not much higher but possibly lower.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end  
    
    if i==37 
Reference="Perry 1985"; 
Species="Katsuwonus pelamis";
Environment="Marine";
Fluid="High-Ventral aorta";
Measurement="Calculated";
Activity="Resting"; 
T_rep=25; 
P_rep=4/760*1e6; %torr to uatm, calculated
P_w_rep=800; %not reported, but "rapidly flowing" through large tank, via seawater intake.
%Assume about same as other tuna experiments at same facility, guess 800
%uatm. i.e., seawater line at this facility seems to have substantial
%respiration before holding tank at the facility's flow rate.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end    
    
    if i==38 
Reference="Perry 1985"; 
Species="Katsuwonus pelamis";
Environment="Marine";
Fluid="High-Ventral aorta";
Measurement="Calculated";
Activity="Active"; 
T_rep=25; 
P_rep=8.7/760*1e6; %torr to uatm, calculated
P_w_rep=800; %not reported, but "rapidly flowing" through large tank, via seawater intake.
%Assume about same as other tuna experiments at same facility, guess 800
%uatm. i.e., seawater line at this facility seems to have substantial
%respiration before holding tank at the facility's flow rate.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end    
    
    if i==39 
Reference="Gonzalez 2001"; 
Species="Amia calva"; %can air-breathe, so using only experiments with no head space and default pH
Environment="Fresh";
Fluid="Low-Dorsal aorta";
Measurement="Calculated";
Activity="Resting"; 
T_rep=15; 
P_rep=0.325/101.325*1e6; %kPa to uatm, calculated
P_w_rep=800; %not reported, but recirculated tap water open to air without 
%depressed water pH, some some unreported aeration scheme likely
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end      
    
    if i==40 
Reference="Gonzalez 2001"; 
Species="Amia calva"; %can air-breathe, so using only experiments with no head space and default pH
Environment="Fresh";
Fluid="Low-Dorsal aorta";
Measurement="Calculated";
Activity="Active"; 
T_rep=15; 
P_rep=0.495/101.325*1e6; %kPa to uatm, calculated
P_w_rep=800; %not reported, but recirculated tap water open to air without 
%depressed water pH, some some unreported aeration scheme likely
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end  
    
    if i==41 
Reference="Turner and Wood 1983"; 
Species="Hippoglossoides elassodon"; 
Environment="Marine";
Fluid="Mod-Caudal artery or vein";
Measurement="Calculated";
Activity="Resting"; 
T_rep=11.5; 
P_rep=2.2/760*1e6; %torr to uatm, calculated
P_w_rep=800; %not reported, pO2 severely undersaturated (0.118 atm!) and  local
%Bamfield BC water typically very high pCO2. Assume about the same as other
%PNW fjord waters, but may be even higher given very low O2
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end     
    
    if i==42 
Reference="Turner and Wood 1983"; 
Species="Hippoglossoides elassodon"; 
Environment="Marine";
Fluid="Mod-Caudal artery or vein";
Measurement="Calculated";
Activity="Active"; %"Wxhausting excercise
T_rep=11.5; 
P_rep=5.3/760*1e6; %torr to uatm, calculated
P_w_rep=800; %not reported, pO2 severely undersaturated (0.118 atm!) and  local
%Bamfield BC water typically very high pCO2. Assume about the same as other
%PNW fjord waters, but may be even higher given very low O2
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end     
    
    if i==43 
Reference="Roth and Rotabakk 2012"; 
Species="Pollachius virens"; 
Environment="Marine";
Fluid="Mod-Caudal vein";
Measurement="Measured";
Activity="Resting"; %caught on line quickly
T_rep=11; 
P_rep=4.4/760*1e6; %mmHg to uatm, calculated
P_w_rep=400; %not reported, but seasonal North Sea seawater ~2008 not far from 400 
%and all excercise in situ
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end    
    
    if i==44 
Reference="Roth and Rotabakk 2012"; 
Species="Pollachius virens"; 
Environment="Marine";
Fluid="Mod-Caudal vein";
Measurement="Measured";
Activity="Active"; %"Exhausting excercise
T_rep=11; 
P_rep=10.8/760*1e6; %mmHg to uatm, calculated
P_w_rep=400; %not reported, but seasonal North Sea seawater ~2008 not far from 400
%and all excercise in situ
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end    
    
    if i==45 
Reference="Schwalme and Mackay 1985"; 
Species="Esox lucius"; 
Environment="Fresh";
Fluid="High-Ventral aorta";
Measurement="Measured";
Activity="Resting"; 
T_rep=19; 
P_rep=5.4/760*1e6; %mmHg to uatm, calculated
P_w_rep=2000; %not reported, but "continuously pumped fresh water". Context
%(Alberta lab facility) suggests tap water with no noted aeration scheme.
%Post hoc, anything less than 1500 uatm would imply unusually high Pmet/Pw ratios,
%which indirectly supports a typical tap value here.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end      

    if i==46 
Reference="Schwalme and Mackay 1985"; 
Species="Esox lucius"; 
Environment="Fresh";
Fluid="High-Ventral aorta";
Measurement="Measured";
Activity="Active"; 
T_rep=19; 
P_rep=11.3/760*1e6; %mmHg to uatm, calculated
P_w_rep=2000; %not reported, but "continuously pumped fresh water". Context
%(Alberta lab facility) suggests tap water with no noted aeration scheme.
%Post hoc, anything less than 1500 uatm would imply unusually high Pmet/Pw ratios,
%which indirectly supports a typical tap value here.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end   

    if i==47 
Reference="Harter 2014"; 
Species="Pterygoplichthys pardalis"; 
Environment="Fresh";
Fluid="Mod-Caudal vein";
Measurement="Calculated";
Activity="Resting"; %"exhaustive excercise, Series 2", note though that
%"excercised" sample was 2 hours after excercise, not immediately
%post-activity as with most studies compiled here
T_rep=28; 
P_rep=(7.2-2.5)/760*1e6; %mmHg to uatm, calculated, no air-breathing observed
P_w_rep=2000; %not reported, but "air-saturated water", freshwater with very 
%low reported pH (6.9); typical Amazon-region fresh waters can be
%substantially higher than 2000uatm, so I view this as a lower bound
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end       
  
    if i==48 
Reference="Harter 2014"; 
Species="Pterygoplichthys pardalis"; 
Environment="Fresh";
Fluid="Mod-Caudal vein";
Measurement="Calculated";
Activity="Active"; %"exhaustive excercise, Series 2", note though that
%"excercised" sample was 2 hours after excercise, not immediately
%post-activity as with most studies compiled here
T_rep=28; 
P_rep=7.2/760*1e6; %mmHg to uatm, calculated, no air-breathing observed
P_w_rep=2000; %not reported, but "air-saturated water", freshwater with very 
%low reported pH (6.9); typical Amazon-region fresh waters can be
%substantially higher than 2000uatm, so I view this as a lower bound
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end    

    if i==49 
Reference="Nelson 1996"; 
Species="Gadus morhua"; 
Environment="Marine"; %Using "full salinity" acclimated Scotian Shelf cod here, at native salinity
Fluid="High-Ventral aorta"; %actually afferant branchial artery, but qualitatively
%similar wrt to expected high vs low CO2.
Measurement="Calculated";
Activity="Resting"; %prior to swimming
T_rep=2; 
P_rep=0.43/101.325*1e6; %kPa to uatm, calculated, no air-breathing observed
P_w_rep=800; %not reported, "flowing filtered seawater" in later winter, Halifax. 
%Note this region has some of the highest seasonal differences in pCO2 in
%the world, so this number is quite uncertain and can only be ballparked
%from other sources (e.g. ship of opportunity transects, such as Becker et
%al. 2018). Unclear how much respiration in seawater line pre-tanks, but
%treating same as other not actively aerated seawater line expts, ~2x local
%water. Post hoc, these Atlantic cod have unusualy high ratios unless
%Pw>1500 uatm, so I view 800 uatm as a lower bound.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end      
    
    if i==50 
Reference="Nelson 1996"; 
Species="Gadus morhua"; 
Environment="Marine"; %Using "full salinity" acclimated Scotian Shelf cod here, at native salinity
Fluid="High-Ventral aorta"; %actually afferant branchial artery, but qualitatively
%similar wrt to expected high vs low CO2.
Measurement="Calculated";
Activity="Active"; %0 h post exhaustion (~same as critical swimming speed)
T_rep=2; 
P_rep=1.07/101.325*1e6; %kPa to uatm, calculated, no air-breathing observed
P_w_rep=800; %not reported, "flowing filtered seawater" in later winter, Halifax. 
%Note this region has some of the highest seasonal differences in pCO2 in
%the world, so this number is quite uncertain and can only be ballparked
%from other sources (e.g. ship of opportunity transects, such as Becker et
%al. 2018). Unclear how much respiration in seawater line pre-tanks, but
%treating same as other not actively aerated seawater line expts, ~2x local
%water. Post hoc, these Atlantic cod have unusualy high ratios unless
%Pw>1500 uatm, so I view 800 uatm as a lower bound.
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end       
    
    if i==51 
Reference="Tang 1992"; 
Species="Oncorhynchus mykiss"; %even more rainbow trout
Environment="Fresh"; 
Fluid="Low-Dorsal aorta"; %Using the measured arterial value (Table 2),
%not the calculate venous values derived from another paper
Measurement="Calculated";
Activity="Resting"; 
T_rep=12; 
P_rep=0.35/101.325*1e6; %kPa to uatm, calculated, no air-breathing observed
P_w_rep=2000; %not reported, more Vancouver tap water 
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end   
    
    if i==52 
Reference="Tang 1992"; 
Species="Oncorhynchus mykiss"; %even more rainbow trout
Environment="Fresh"; 
Fluid="Low-Dorsal aorta"; %Using the measured arterial value (Table 2),
%not the calculate venous values derived from another paper
Measurement="Calculated";
Activity="Active"; %excercised to exhaustion (chased in tank)
T_rep=12; 
P_rep=1.02/101.325*1e6; %kPa to uatm, calculated, no air-breathing observed
P_w_rep=2000; %not reported, more Vancouver tap water 
%Report out
Temperature_Reported=T_rep;pCO2_Reported=P_rep;pCO2_water=P_w_rep;
%Add to table
Pint(i,:)=table(Reference,Species,Environment,Fluid,Activity,Measurement,Temperature_Reported,...
    pCO2_Reported,pCO2_water);
clearvars T_rep P_rep P_w_rep;
clearvars Reference Species Fluid Measurement Temperature_Reported Activity Environment;
clearvars pCO2_Reported pCO2_water;
    end     
    
end  
save BloodPint.mat Pint -v7.3
else
    load BloodPint.mat
end


%%

A=Pint{:,[8,9]};
B=(A(:,1)-A(:,2))./A(:,2);
%All
subplot(221);histogram(B(:,1),[0:1:30]);xlabel('Pmet/Pw');ylabel('Meas');
legend('All');
%Active+Resting paired only
il=[8,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51]; %note 7 high, 8 low, for this pair only
ih=[7,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52];%don't pair 5,6, Pw too different
iil=zeros(length(A),1);iil(il)=1;iil=logical(iil);
iih=zeros(length(A),1);iih(ih)=1;iih=logical(iih);
subplot(223);histogram(B(iil,1),[0:1:30]);hold on;histogram(B(iih,1),[0:1:30]);xlabel('Pmet/Pw');ylabel('Meas');
legend('Paired resting','Paired active');
subplot(224);histogram(B(iih,1)./B(iil,1),[0:0.5:4]);xlabel('Active/Resting');ylabel('Meas');
%Circulation: expected low, moderate, or high pCO2
ibh=[2,19,20,31,32,45,46,49,50];
iibh=zeros(length(A),1);iibh(ibh)=1;iibh=logical(iibh);
subplot(222);histogram(B(~iibh,1),[0:1:30]);hold on;histogram(B(iibh,1),[0:1:30]);
xlabel('Pmet/Pw');ylabel('Meas');
legend('All arterial','All venous');

%Expected whole blood mean after assumptions about blood source
