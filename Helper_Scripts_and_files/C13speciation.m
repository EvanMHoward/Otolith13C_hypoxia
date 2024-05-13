function [R13CO2g R13CO2aq R13HCO3 R13CO3 R13Arag R13DIC]=...
        C13speciation(R13known,pH,K1,K2,T,speciesknown)

%Author: Evan Howard (ehoward2@uw.edu)
%Original version: December 2020 (modifed from 'd13CbySpecies.m' by Alex
%Gagnon (University of Washington School of Oceanography)
%Current version: Februrary 2022

%Description: %This recalculates equilibrium constants and associated 
% equilibrium isotopic fractionations associated with that speciation for 
% alternative pH (total scale) and ionic strength solutions. See notes in
% 'carbeq.m'.


%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 

%The ratio R=(13C/12C)_sample (absolute, already multipled by standard)
%of a particular species, pH and K values on total pH scale (NOT SEAWATER 
%OR NBS), T (degrees Celsius), S (PSS-78),
%and the carbon species that has been measured.

%What is already known?
if strcmp(speciesknown,'DIC')==1 %In all cases these should already be absolute ratios
    R13_DIC=R13known;qq=R13_DIC;
elseif strcmp(speciesknown,'CO2g')==1
    R13_CO2g=R13known;qq=R13_CO2g;
elseif strcmp(speciesknown,'CO2aq')==1
    R13_CO2aq=R13known;qq=R13_CO2aq;
elseif strcmp(speciesknown,'HCO3')==1
    R13_HCO3=R13known;qq=R13_HCO3;
elseif strcmp(speciesknown,'CO3')==1
    R13_CO3=R13known;qq=R13_CO3;
elseif strcmp(speciesknown,'Arag')==1
    R13_Arag=R13known;qq=R13_Arag;
end

%What are mole fractions of each species in DIC?
H=10.^(-1.*pH);D=(H.^2)+H.*K1+K1.*K2;
X0=(H.^2)./D;X1=H.*K1./D;
X2=K1.*K2./D;
%X2 is the term most likely to in error in a solution of very different 
% ionic strength, so instead could solve by difference (X2=1-X0-X1). 
% However, if sodium carbonate salts are appreciable (usually at high pH) 
% then this substantially overestimates the CO3 component.

%What are fractionation factors between species?
%From Zhang Quay, & Wilbur 1995 (https://doi.org/htt10.1016/0016-7037(95)91550-D)
%NOTE: The relationship presented for DIC-gas is valid over ph=7.5 to 8.5 
% (see table in publication), %which encompasses most fish and ocean 
% applications (actual validity is related to fraction of CO3, 0.05-0.2 
% range, which spans a larger pH range in fish). However, these ranges and 
% assumptions are not required if estimate R_DIC directly from species 
% fractionation factors.

%isotope effects
 E_HCO3_CO2g=-(0.1141).*T+(10.78); %Epsilon, bicarbonate produced from CO2 gas
 E_CO2aq_CO2g=0.0049.*T-1.31; %Epsilon, CO2 aq produced from CO2 gas
 E_CO3_CO2g=-0.052.*T+7.22; %Epsilon, CO3 produced from CO2 gas
E_Arag_HCO3=2.7; %Romanek,Grossman,&Morse (https://doi.org/10.1016/0016-7037(92)90142-6)

%fractionation factors: as with all fractionation calculations, using 
% fractionation factors (alphas) limits numerical errors introduced when
% adding isotope effects (epsilons)
a_HCO3_CO2g=E_HCO3_CO2g./1000+1;
a_CO2aq_CO2g=E_CO2aq_CO2g./1000+1;
a_CO3_CO2g=E_CO3_CO2g./1000+1;
a_Arag_HCO3=E_Arag_HCO3./1000+1;

if strcmp(speciesknown,'DIC')==1
    %Zero-finder only for the case when R_DIC is known
    mbslv= @(x)R13massbalsolver(x,X0,X1,X2,a_HCO3_CO2g,a_CO2aq_CO2g,a_CO3_CO2g);
    R13_HCO3= fzero(mbslv,qq);
    R13_CO2g=R13_HCO3.*(a_HCO3_CO2g.^-1);
    R13_CO2aq=R13_HCO3.*(a_HCO3_CO2g.^-1).*(a_CO2aq_CO2g);
    R13_CO3=R13_HCO3.*(a_HCO3_CO2g.^-1).*(a_CO3_CO2g);
    R13_Arag=R13_HCO3.*a_Arag_HCO3;
elseif strcmp(speciesknown,'CO2g')==1
    R13_HCO3=qq.*(a_HCO3_CO2g);
    R13_CO2aq=qq.*(a_CO2aq_CO2g);
    R13_CO3=qq.*(a_CO3_CO2g);
    R13_Arag=R13_HCO3.*a_Arag_HCO3;  
    R13_DIC=X0.*R13_CO2aq+X1.*R13_HCO3+X2.*R13_CO3;  
elseif strcmp(speciesknown,'CO2aq')==1
    R13_CO2g=qq.*(a_CO2aq_CO2g.^-1);
    R13_HCO3=R13_CO2g.*(a_HCO3_CO2g);
    R13_CO3=R13_CO2g.*(a_CO3_CO2g);
    R13_Arag=R13_HCO3.*a_Arag_HCO3;  
    R13_DIC=X0.*R13_CO2aq+X1.*R13_HCO3+X2.*R13_CO3;
elseif strcmp(speciesknown,'HCO3')==1
    R13_CO2g=qq.*(a_HCO3_CO2g.^-1);
    R13_CO2aq=R13_CO2g.*(a_CO2aq_CO2g);
    R13_CO3=R13_CO2g.*(a_CO3_CO2g);
    R13_Arag=R13_HCO3.*a_Arag_HCO3;  
    R13_DIC=X0.*R13_CO2aq+X1.*R13_HCO3+X2.*R13_CO3;    
elseif strcmp(speciesknown,'CO3')==1
    R13_HCO3=qq.*(a_CO3_CO2g.^-1).*(a_HCO3_CO2g);
    R13_CO2g=R13_HCO3.*(a_HCO3_CO2g.^-1);
    R13_CO2aq=R13_CO2g.*(a_CO2aq_CO2g);
    R13_Arag=R13_HCO3.*a_Arag_HCO3;  
    R13_DIC=X0.*R13_CO2aq+X1.*R13_HCO3+X2.*R13_CO3;  
elseif strcmp(speciesknown,'Arag')==1
    R13_HCO3=qq.*(a_Arag_HCO3.^-1);
    R13_CO2g=R13_HCO3.*(a_HCO3_CO2g.^-1);
    R13_CO2aq=R13_CO2g.*(a_CO2aq_CO2g);
    R13_CO3=R13_CO2g.*(a_CO3_CO2g);
    R13_DIC=X0.*R13_CO2aq+X1.*R13_HCO3+X2.*R13_CO3;     
end
%OUT=[R13_CO2g,R13_CO2aq,R13_HCO3,R13_CO3,R13_Arag,R13_DIC];
R13CO2g=R13_CO2g;R13CO2aq=R13_CO2aq;R13HCO3=R13_HCO3;
R13CO3=R13_CO3; R13Arag=R13_Arag; R13DIC=R13_DIC;

%% Zero-finder function (from Alex Gagnon)
function errval = R13massbalsolver(R13_HCO3,X0,X1,X2,a_HCO3_CO2g,a_CO2aq_CO2g,a_CO3_CO2g)
    R13_CO2g=R13_HCO3.*(a_HCO3_CO2g.^-1);
    R13_CO2aq=R13_HCO3.*(a_HCO3_CO2g.^-1).*(a_CO2aq_CO2g);
    R13_CO3=R13_HCO3.*(a_HCO3_CO2g.^-1).*(a_CO3_CO2g);
    R13_Arag=R13_HCO3.*a_Arag_HCO3;
%calculate abundances
A_DIC=R13_DIC/(R13_DIC+1);
A_CO2aq=R13_CO2aq/(R13_CO2aq+1);
A_HCO3=R13_HCO3/(R13_HCO3+1);
A_CO3=R13_CO3/(R13_CO3+1);
%calcuate mass balance zero function
errval = A_CO2aq*X0+A_HCO3.*X1+A_CO3.*X2-A_DIC;
end

end
