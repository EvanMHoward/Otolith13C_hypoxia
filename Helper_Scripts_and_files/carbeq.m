function [K0 K1 K2]=carbeq(T,S)

%Author: Evan Howard (ehoward2@uw.edu)
%Original version: December 2020
%Current version: September 2022

%Description: This function calculates carbonate system equilibria. 
%Carbonate system equilibria K1 and K2 are derived for use with the pH_Tot
%scale using Millero 2010 (doi:10.1071/MF09254). K0 (Henry's law
%coefficient is derived using the equations of Weiss 1974, as described in
%Dickson et. al 2007 (Dickson, Sabine, Christian, 2007, Guide to Best
%Practices for Ocean CO2 Measurements. PISCES Special Publication 3). 

%The total pH scale differs from the NBS scale typically used by 
% physiologists or the the seawater pH scale previosuly in widespread use
% in oceanography. The difference between the total scale and the former 
% is ~-0.15, while that with the latter is ~+0.01. This can be exactly
% computed using CO2SYS, but for the purposes of this work the
% approximate differences are sufficient. Teleost fish plasma and endolymph
% have ~1/3 the ionic strength of seawater. When using this script an 
% isoionic adjustments can be applied to input salinity, but a more robust 
% approach is to revise the activity coefficients using the Davies 
% formulation of the Debye-Huckel theory.

%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 

Tabs=T+273.15;

lnK0=93.4517.*(100./Tabs)-60.2409+23.3585.*log(Tabs./100)+...
    S.*(0.023517-0.023656.*(Tabs./100.)+0.0047036.*((Tabs./100).^2));
K0=exp(lnK0); %where K0=[CO2(aq)]/f(CO2) (mol kg-1 atm-1 of solution)

%NOTE: this gives nearly identical results to CO2SYS but is much faster, 
% using total pH scale which is slightly different in K2 from 
% seawater scale but not coded for Millero 2010 option in CO2SYS.

%e.g. CO2SYSV3(SWS)->0.0371001120459998	9.10947386760744e-07 4.20308124035023e-10
%and this script->   0.0371001120459998	9.09764638729761e-07 4.16271134879162e-10

pK1n=-126.34048+6320.813./Tabs+19.568224.*log(Tabs);
pK2n=-90.18333+5143.692./Tabs+14.613358.*log(Tabs);
a1=[13.4051;0.03185;-5.218e-5;-531.095;-5.7789;-2.0663];
a2=[21.5724;0.1212;-3.714e-4;-798.292;-18.951;-3.403];

a=a1;pkn=pK1n;
pK=(a(1).*(S.^0.5)+a(2).*S+a(3).*(S.^2))+...
    (a(4).*(S.^0.5)+a(5).*S)./Tabs+...
    (a(6).*(S.^0.5)).*log(Tabs)+pkn;
K1=10.^(-pK);
a=a2;pkn=pK2n;
pK=(a(1).*(S.^0.5)+a(2).*S+a(3).*(S.^2))+...
    (a(4).*(S.^0.5)+a(5).*S)./Tabs+...
    (a(6).*(S.^0.5)).*log(Tabs)+pkn;
K2=10.^(-pK);
end