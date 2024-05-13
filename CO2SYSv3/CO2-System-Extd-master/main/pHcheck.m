PAR1=2053;
PAR2=8.00;
PAR1TYPE=1;
PAR2TYPE=3;
SAL=30;
TEMPIN=20;TEMPOUT=TEMPIN;
PRESIN=0;PRESOUT=PRESIN;
SI=20;PO4=1;NH4=1;H2S=0; %Puget Sound approximate values
pHSCALEIN=4;
K1K2CONSTANTS=15; %Seawater scale, but almost same as total and more appropriate to natural waters
KSO4CONSTANT=3; %for consistency with above
KFCONSTANT=2;
BORON=2;

A=CO2SYS(PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL, ...
    TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,NH4,H2S,...
    pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON);

%    39 - pH input (Total)     ()          
%    40 - pH input (SWS)       ()          
%    41 - pH input (Free)      ()          
%    42 - pH input (NBS)       () 

PAR2 %This went in
[A(39),A(40),A(42)] %This comes out, Tot, SWS, NBS

%Now change Alk and compare
PAR1b=PAR1+127.4;
B=CO2SYS(PAR1b,PAR2,PAR1TYPE,PAR2TYPE,SAL, ...
    TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,NH4,H2S,...
    pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON);

[B(39),B(40),B(42)] %This comes out, Tot, SWS, NBS

%Now change S and compare
SALb=25;
C=CO2SYS(PAR1,PAR2,PAR1TYPE,PAR2TYPE,SALb, ...
    TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,NH4,H2S,...
    pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON);

[C(39),C(40),C(42)] %This comes out, Tot, SWS, NBS