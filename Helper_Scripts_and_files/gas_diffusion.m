function [diff_Ne,diff_Ar,diff_Kr,diff_Xe,diff_O2,diff_N2,diff_He] = gas_diffusion(S,T)

% gas_diffusion   Diffusion coefficients for various gases in fresh/sea water
%=========================================================================
% gas_diffusion Version 2.0 16 July 2013
%          Author: Roberta C. Hamme (University of Victoria)
%
% USAGE:  [diff_Ne,diff_Ar,diff_Kr,diff_Xe,diff_O2,diff_N2,diff_He] = gas_diffusion(S,T)
%
% DESCRIPTION:
%    Diffusion coefficients of various gases in fresh/sea water
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = salinity    [PSS-78]
%   T = temperature [degree C]
%
% OUTPUT:
%   diff_Ne = diffusion coefficient of Ne  [cm^2/s or 10^-4 m^2/s] 
%   diff_Ar = diffusion coefficient of Ar  [cm^2/s or 10^-4 m^2/s] 
%   diff_Kr = diffusion coefficient of Kr  [cm^2/s or 10^-4 m^2/s] 
%   diff_Xe = diffusion coefficient of Xe  [cm^2/s or 10^-4 m^2/s] 
%   diff_O2 = diffusion coefficient of O2  [cm^2/s or 10^-4 m^2/s] 
%   diff_N2 = diffusion coefficient of N2  [cm^2/s or 10^-4 m^2/s] 
%   diff_He = diffusion coefficient of He  [cm^2/s or 10^-4 m^2/s] 
% 
% AUTHOR:  Roberta Hamme (rhamme@uvic.ca)
%
% REFERENCE:
%    He, Ne, Kr, Xe freshwater values from Jahne et al., 1987.
%       "Measurement of Diffusion Coeffients of Sparingly Soluble Gases in Water"
%       J. Geophys. Res., 92(C10), 10767-10776.
%    Ar freshwaters values are extrapolated from Jahne et al. 1987
%       He, Ne, Kr, Xe values at each temperature were fitted to D vs. mass^-0.5 
%       relationship to predict Ar at those temperatures, then Ar was fit to a 
%       ln(D_Ar) vs. 1/T(K) relationship to obtain Eyring equation coefficients
%    O2 and N2 freshwater values from  Ferrell and Himmelblau, 1967.
%       "Diffusion coefficients of nitrogen and oxygen in water"
%       J. Chem. Eng. Data, 12(1), 111-115, doi: 10.1021/je60032a036.
%    Correction for salinity is based on Jahne's observed average 4.9% decrease in 
%       diffusivity for H2 and He in 35.5 ppt NaCl solution
%
%    for Ne, the Jahne values compare well with and fall between those of
%       Wise and Houghton 1968 and Holz et al. 1994
%    for Ar, the extrapolated Jahne values compare well with Wise and Houghton 1968,
%       O'Brien and Hyslop 1977, and a numerical simulation by Bourg et al. 2008
%       but are higher than other reported values
%    for Kr, the Jahne values compare well with Wise and Houghton 1968,
%       and a numerical simulation by Bourg et al. 2008
%    for Xe, the Jahne values compare well with Pollack 1981, and a numerical 
%       simulation by Bourg et al. 2008, but fall significantly above Wise and Houghton 1968
%       and below Weingartner et al. 1992
%    for O2, there is general agreement among measurements. The Ferrel and Himmelblau values 
%       agree reasonably well with Baird and Davidson 1962, Wise and Houghton 1966,
%       Duda and Vrentas 1968, O'Brien and Hyslop 1977, and the Wilke and Change (1955) theory 
%       as tabulated by Wanninkhof 1992, but lie below Krieger et al 1967
%    for N2, there is less agreement. The Ferrel and Himmelblau values 
%       agree reasonably well with Baird and Davidson 1962, O'Brien and Hyslop 1977, 
%       and the Wilke and Change (1955) theory as tabulated by Wanninkhof 1992, 
%       but lie significantly below the values of Wise and Houghton 1966 and Krieger et al 1967
%    for He, I did not investigate comparisons of data, but chose Jahne 
%       since their work for other gases appears to be the best
%
% DISCLAIMER:
%    This software is provided "as is" without warranty of any kind.  
%=========================================================================

% CALLER: general purpose
% CALLEE: none


%----------------------
% Check input parameters
%----------------------
if nargin ~=2
   error('gas_diffusion.m: Must pass 2 parameters')
end %if

% CHECK S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);

  
% Check that T&S have the same shape or are singular
if ((ms~=mt) | (ns~=nt)) & (ms+ns>2) & (mt+nt>2)
   error('gas_diffusion.m: S & T must have same dimensions or be singular')
end %if

%------
% BEGIN
%------

% Calculate temperature in Kelvin
TK = T + 273.15;

% Apply Eyring equation and appropriate coefficients for each gas in freshwater
diff_Ne_0sal = 1608 * exp(-14.84./(8.3144621e-3 * TK)) * 1e-5;         
diff_Ar_0sal = 2227 * exp(-16.68./(8.3144621e-3 * TK)) * 1e-5;         
diff_Kr_0sal = 6393 * exp(-20.20./(8.3144621e-3 * TK)) * 1e-5;         
diff_Xe_0sal = 9007 * exp(-21.61./(8.3144621e-3 * TK)) * 1e-5;         
diff_N2_0sal = 3412 * exp(-18.50./(8.3144621e-3 * TK)) * 1e-5;         
diff_O2_0sal = 4286 * exp(-18.70./(8.3144621e-3 * TK)) * 1e-5;         
diff_He_0sal =  818 * exp(-11.70./(8.3144621e-3 * TK)) * 1e-5;         

% Apply salinity correction of 4.9% per 35.5 ppt salinity
diff_Ne = diff_Ne_0sal .* (1 - 0.049 * S / 35.5);
diff_Ar = diff_Ar_0sal .* (1 - 0.049 * S / 35.5);
diff_Kr = diff_Kr_0sal .* (1 - 0.049 * S / 35.5);
diff_Xe = diff_Xe_0sal .* (1 - 0.049 * S / 35.5);
diff_O2 = diff_O2_0sal .* (1 - 0.049 * S / 35.5);
diff_N2 = diff_N2_0sal .* (1 - 0.049 * S / 35.5);
diff_He = diff_He_0sal .* (1 - 0.049 * S / 35.5);

return