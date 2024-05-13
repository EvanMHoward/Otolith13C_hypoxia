%estimated values from PRISM bottle data from mid-2000s
%http://www.prism.washington.edu/story/Observing+the+Sound

TA=2000;ePAR1=50; %vector or value of alkalinity to check
pH=[7.5:0.05:8.5];ePAR2=0.008; %vector of pH to check
S=20;eSAL=0; %salinity, and error is bias rather than random so ignored here
Ti=20; Tf=Ti; eTEMP=1; %temperature, initial and final (in lab), %half degree precision uncertainty based on instrument agreement
Pi=0; Pf=Pi; %pressure, initial and final
SI=20; eSI=0; %set these estimated values to 0, any bias is not random error in this case
PO4=1; ePO4=0;
NH4=1;eNH4=0;
H2S=0;eH2S=0;
eCa=0;eB=0;
r=0; %correlation coefficient between input carbonate system parameters, assumed
pHsc=1; %total scale
K1K2=10; %10 is Leuker 2000,  total pH over appropriate ranges. 
% Note that the 2014 waters and millero paper has total pH scale too, just 
% needs to be added as case to this script. Also fine to use Seawater
% scale, which is very close, and actually may be better measured in these
% options and more consistent with other partition coefficients.
kSO4=1; %Dickson, but Waters and Millero may be better particularly if used with SW scale.
kF=2; 
TB=2; %Lee et al 2010


%
%  **** SYNTAX:
%  [err, headers, units] = errors(PAR1,PAR2,PAR1TYPE,PAR2TYPE,...
%                                 SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,...
%                                 NH4,H2S,ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,...
%                                 eNH4,eH2S,epK,eBt,r,pHSCALEIN,K1K2CONSTANTS,...
%                                 KSO4CONSTANT,KFCONSTANT,BORON,eCAL(optional))

H=nan(length(TA),length(pH));
for i=1:length(TA)
    ta=TA(i);
    for j=1:length(pH)
    ph=pH(j);
        [Result]=errors(ta,ph,1,3,S,Ti,Tf,Pi,Pf,SI,PO4,NH4,H2S,...
            ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,eNH4,eH2S,[],eB,r,pHsc,K1K2,...
            kSO4,kF,2,eCa);
        H(i,j)=Result(3);
    end
end

%%
%H is in nmol kg-1, so multiply by 10^-9 to get H+ concentration in moles
%for comparison of pH error bounds.
%Example where TA=2000 umolal
pHbnd=nan(3,length(pH));
pHbnd(1,:)=(-log10((10.^(-pH))-H(1,:).*(10^-9))); %change row of H if more than one TA tested
pHbnd(2,:)=pH;
pHbnd(3,:)=(-log10((10.^(-pH))+H(1,:).*(10^-9)));

figure;
plot(pHbnd(2,:),'-k');hold on;
plot(pHbnd(1,:),'-r');
plot(pHbnd(3,:),'-r');

pHunc=-log10((10.^(-pH))-H.*(10^-9))-pH;open pHunc