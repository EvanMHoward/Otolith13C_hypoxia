function [Dco2,Dhco3,Dco3]=co2_diffusion(S,T)

%Author: Evan Howard (ehoward2@uw.edu)
%Original version: December 2020
%Current version: November 2021

%Description: %This function calculate the diiffusion coefficients of
%carbonate system chemical species. This is based on the freshwater
% parameterization of Zeebe 2011 (doi: 10.1016/j.gca.2011.02.010), scaled 
% with the same salinity dependence as Roberta Hamme's gas_diffusion.m
%script. This may introcuce errors due to the much higher solubility of
%CO2, carbonate system contribution to salinity, and different ionic
%interactions for dissolved carbonate system species, but is the most
%straightforward approach given the lack of experimental parameterizations
%for seawater.

%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 

%units, 10-9 m2 s-1
Dco2_0sal = 14.6836.*(((T+273.15)./217.2056-1).^1.9970);
Dhco3_0sal = 7.0158.*(((T+273.15)./204.0282-1).^2.3942);
Dco3_0sal = 5.4468.*(((T+273.15)./210.2646-1).^2.1929);

%convert to cm2 s-1, or 10-4 m2 s-1
Dco2_0sal=Dco2_0sal.*(1e-9)./(1e-4);
Dhco3_0sal=Dhco3_0sal.*(1e-9)./(1e-4);
Dco3_0sal=Dco3_0sal.*(1e-9)./(1e-4);

%scale for salinity versus freshwater: 
%This is appropriate for an approximately spherical ideal gas (Argon), but 
% has errors for these nonideal solutes that are beyond the scope of this 
%script. Nonetheless, it seems to do a very good job matching the offset 
% between experimental seawater and freshwater data (or freshwater molecular
%dynamics simulations) for CO2 and HCO3 (no data for CO3) based on
%literature review.
Dco2 = Dco2_0sal.*(1-0.049.*S./35.5)';
Dhco3 = Dhco3_0sal.*(1-0.049.*S./35.5)'; 
Dco3 = Dco3_0sal.*(1-0.049.*S./35.5)';
end