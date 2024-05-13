function [pH]=pH_T(T,S,pHref,Tref)

%Author: Evan Howard (ehoward2@uw.edu)
%Original version: December 2020
%Current version: December 2020

%Description: %pH changes as a function of temperature. Reference and 
% output must be on total pH scale (~seawater scale). 

%Uses the temperature dependence of Kw (Millero 1995, as in Dickson et al.
%2007 (Dickson, Sabine, Christian, 2007, Guide to Best
%Practices for Ocean CO2 Measurements. PISCES Special Publication 3), to 
%recalculate pH for the same solution with warming or cooling. 
%The temperature dependence of pH is observed to parallel that of Kw with 
%a constant relative alkalinity (unchanged acidity because of concurrent 
%changes in bases) in both the ocean 
%(Hunter 1998, https://doi.org/10.1016/S0967-0637(98)00047-8) and
%fish blood (Cameron 1998, https://doi.org/10.1016/0034-5687(78)90092-0),
%as expected from pure thermodynamic considerations. In practice, this
%translates to a -0.01 to -0.02 change in pH for each 1 degree increase in
%temperature.

%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 

%temperature dependence
Tabs=T+273.15;
lnKw=(-13847.26./Tabs+148.9652-23.6521.*log(Tabs)+...
    (118.67./Tabs-5.977+1.0495.*log(Tabs)).*(S.^0.5)-0.01615.*S);
Tabs=Tref+273.15;
lnKwref=(-13847.26./Tabs+148.9652-23.6521.*log(Tabs)+...
    (118.67./Tabs-5.977+1.0495.*log(Tabs)).*(S.^0.5)-0.01615.*S);

%change scales
pKw=-log10(exp(lnKw));pKwref=-log10(exp(lnKwref));

%calculate new pH
pH=(pKw-pKwref)+pHref;
end