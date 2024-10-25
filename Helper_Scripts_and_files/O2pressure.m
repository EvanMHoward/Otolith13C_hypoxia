function po2 = O2pressure(O2,T,S,z)
% O2PRESSURE This function compute partical pressure of oxygen.
%
% Author: unknown
% Original version: unknown
% Current version: August 27, 2024

%Description: Hydrostatic correction for the increasing pO2 of seawater
%with depth. Sources and code performance checked.

%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 


%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes partial pressure of O2 in seawater based on 
% Inputs: O2 concentration (umol/kg), T, S, and depth
% Includes correction for the effect of hydrostatic pressure of the water column.
% based on Enns et al., J. Phys. Chem. 1964
%  d(ln p)/dP = V/RT
% where p = partial pressure of O2, P = hydrostatic pressure
% V = partial molar volume of O2, R = gas constant, T = temperature
%

%EH addendum:
%The solubility of O2 (though not other gases) increases <0.1% per 1000 m
%depth (Klotz 1963; doi:10.4319/lo.1963.8.2.0149). This can be neglected 
%compared to the much larger effect of increasing pressure on the partial 
%pressure in equilibrium with any given concentration (Enns et al. 1964; 
%doi:10.1021/j100886a005) See Eckert (1973; 10.1126/science.180.4084.426)
%for additional discussion of how the same phenonmenon results in both
%effects.  

%%%%%%%%%%%%%%%%%%%%%%

%% constants
XiO2=0.209; % Mean atmospheric O2 mixing ratio
Patm=1; % Atm pressure
V=32e-6; % partial molar volume of O2 (m3/mol)
R=8.31; % Gas constant [J/mol/K]
db2Pa=1e4; % convert pressure: decibar to Pascal

%% Solubility with pressure effect

P=sw_pres(z,z*0); % seawater pressure [db] !! Warning - z*0 neglects gravity differences w/ latitude

dP=P*db2Pa;
pCor=exp(V*dP./(R*(T+273.15))); % 

Kh=O2sol(S,T)/(Patm*XiO2); % solubility [umol/kg/atm]
po2=(O2./Kh).*pCor;






