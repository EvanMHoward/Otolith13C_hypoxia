%Author: Evan Howard (ehoward2@uw.edu)
%Original version: January 2021
%Current version: March 2021

%Description: This script extracts isotope and metadata for Pacific cod
%from Helser et al. 2014 database.

%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 

%% Import isotopic data
%Didn't automate this, manually load data from Isotopes.xlsx

%Import Isotopes.xlsx, table 2, columns A,B,E-G,M,O 
%Name this resulting mat field IsotopesS1

Isotopes=IsotopesS1;Isotopes(1,:)=[];
Isotopes2=Isotopes;

%DO NOT remove missing data, still need the width information
ID=Isotopes2{:,2};
Sample=Isotopes2{:,1};

Haul=nan(size(ID,1),1);Cruise=nan(size(ID,1),1);
for i=1:length(Haul)
    tmp=num2str(ID(i));
    tmp2s=char(extractBetween(tmp,1,3));
    tmp2=str2double(tmp2s);
    if str2double(tmp2s(1))==1
        Haul(i)=tmp2;
    else
        Haul(i)=str2double(tmp2s(1:2));
    end
    if Haul(i)<100
       Cruise(i)=str2double(extractBetween(tmp,3,8));
    else
       Cruise(i)=str2double(extractBetween(tmp,4,9));
    end
end

%Manually correct missing ages caused by inconcsistent notation if needed
Age=nan(size(ID,1),1);
for i=1:length(ID)
tmp=Isotopes2{i,3};
if strcmp(tmp,"Year 1")|strcmp(tmp,"year 1")|strcmp(tmp,"Mark 1")|strcmp(tmp,"Year 1 start new 1"); Age(i)=1;
elseif strcmp(tmp,"Year 2")|strcmp(tmp,"year 2")|strcmp(tmp,"Year 2?"); Age(i)=2;
elseif strcmp(tmp,"Year 3")|strcmp(tmp,"year 3")|strcmp(tmp,"Year 3-Check"); Age(i)=3;
elseif strcmp(tmp,"Year 4")|strcmp(tmp,"year 4"); Age(i)=4;
elseif strcmp(tmp,"Year 5")|strcmp(tmp,"year 5"); Age(i)=5;
elseif strcmp(tmp,"Year 6")|strcmp(tmp,"year 6"); Age(i)=6;
end
end

IsotopeInfo=[Cruise,Haul,Sample,Age,Isotopes2{:,4:7}];

%save('IsotopeInfo','IsotopeInfo'); %uncomment if saving