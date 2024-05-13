%Author: Evan Howard (ehoward2@uw.edu)
%Original version: March 2022
%Current version: April 2024

%Description: This script processes and plots the modeled global otolith
%distributions generated from 'HowardDeutsch_Global_....m'. It also
%produces some associated supplementary figures.

%No express or implied warranty: This script, as well as any supplementary 
% scripts or data tables provided for its operation, is provided "as is"
% without warranty or guarantee of any kind. 

%% Prepare workspace
addpath('.../Helper_Scripts_and_files/'); %Location of helper scripts
load Goc; %WOA gridding
load InterSpecies_b.mat; %Isotope observations
load d3D.mat; %modeled global d13Coto distributions
%% End-member d13C and pCO2 distributions (Fig 1)
set(0,'DefaultAxesFontSize', 12);
set(0,'defaultLineLineWidth', 2);

dEdge=-25:5;dCent=dEdge(2:end)-.5;
pEdge=0:100:12000;pCent=pEdge(2:end)-100;
%Kernals for cdf smoothing (both 10% of plotted range)
kn1=0.1;
kn2=0.1;
stn=0.011180; %VPDB standard

%Observed from small otolith dataset
obsoto=nan(length(C_IS),100);
obsdic=nan(length(C_IS),100);
obsom=nan(length(C_IS),100);
obsco2=nan(length(C_IS),100);
aa=nan(length(C_IS),100);
for i=1:length(C_IS) %13-16 are tropical species with no colocated 
    % organic matter data, estimated tissues from other animals that are in
    % remote geographic areas
    tmp=C_IS{i,:};
    t=tmp(:,2); %use to recalculate dic to CO2aq assuming pH=8.1 and S=34
    o=tmp(:,1);d=tmp(:,4);om=tmp(:,5);
    [K0 K1 K2]=carbeq(t,34.*ones(length(t),1));
    Ref=((d)./1000+1).*stn;
    Rc=nan(length(t),1);
    for j=1:length(t)
        if isnan(Ref(j)); Rc(j,1)=nan;
        else
    [~,Rci,~,~,~,~]=C13speciation(Ref(j),8.1,K1(j),K2(j),t(j),'DIC'); %CO2aq,seawater
    Rc(j,1)=Rci;clearvars Rci
        end
    end
    c=(Rc./stn-1).*1000;
    obsoto(i,1:length(o))=o;
    obsdic(i,1:length(d))=d;
    obsco2(i,1:length(c))=c;
    obsom(i,1:length(om))=om;
    aa(i,1:length(om))=om-c;
    clearvars tmp o d om Rc c Ref K0 K1 K2 j;
end
%Warning, scaling will be wrong if don't excise nan.
%Irrelevant if solving for probability from pdf instead, but if need to
%plot raw histograms be sure to index and call histcounts only for the
%numerical values.
[mN_delOto]=histcounts(obsoto(:),dEdge,'Normalization','probability'); 
[mN_delDIC]=histcounts(obsdic(:),dEdge,'Normalization','probability');
[mN_delCO2]=histcounts(obsco2(:),dEdge,'Normalization','probability');
[mN_delOM]=histcounts(obsom(:),dEdge,'Normalization','probability');
pdfdoto=fitdist(obsoto(:),'Kernel','Kernel','normal','Support','unbounded','Width',kn1);
pdfddic=fitdist(obsdic(:),'Kernel','Kernel','normal','Support','unbounded','Width',kn1);
pdfdco2=fitdist(obsco2(:),'Kernel','Kernel','normal','Support','unbounded','Width',kn1);
pdfdom=fitdist(obsom(:),'Kernel','Kernel','normal','Support','unbounded','Width',kn1);
cdoto=cdf(pdfdoto,dEdge);cddic=cdf(pdfddic,dEdge);cdco2=cdf(pdfdco2,dEdge);cdom=cdf(pdfdom,dEdge);
clearvars pDoto pDdic pDco2 pDom;
for i=1:(length(dEdge)-1) %calculate probabilities associated with every 100 uatm bin based on pdfs
    pDoto(i)=cdoto(i+1)-cdoto(i);
    pDdic(i)=cddic(i+1)-cddic(i);
    pDco2(i)=cdco2(i+1)-cdco2(i);
    pDom(i)=cdom(i+1)-cdom(i);
end

figure(1);
subplot(2,1,1);
plot(dCent,pDom,'Color',[1 (165/255) 0]);hold on;
plot(dCent,pDco2,'b');plot(dCent,pDoto,'k');plot(dCent,dCent.*0,'k');
xlabel('\delta^{13}C (%\fontsize{8}o)'); ylabel('{\it{f}}_{observed}');
ylim([0 0.45]);xlim([-25 5]);
yticks([0:0.1:0.4]);xticks([-20:10:0]);
set(gca,'TickLength',[0.01 0.01],'LineWidth',2,'box','off','Layer','top');

%Observed pCO2 distributions
load BloodPint.mat;
A=Pint{:,[8,9]};
[mN_Pw]=histcounts(A(:,2),pEdge,'Normalization','probability');
[mN_Pint]=histcounts(A(:,1),pEdge,'Normalization','probability');
pdfPw=fitdist(A(:,2),'Kernel','Kernel','normal','Support','positive','Width',kn2);
pdfPint=fitdist(A(:,1),'Kernel','Kernel','normal','Support','positive','Width',kn2);
cPw=cdf(pdfPw,pEdge);cPint=cdf(pdfPint,pEdge);clearvars pPw pPint;
for i=1:(length(pEdge)-1) %calculate probabilities associated with every 100 uatm bin based on pdfs
    pPw(i)=cPw(i+1)-cPw(i);
    pPint(i)=cPint(i+1)-cPint(i);
end
subplot(2,1,2);
plot(pCent,pPw,'-b','LineWidth',2);hold on;plot(pCent,pPint,'-r','LineWidth',2);
plot(pCent,pCent.*0,'k');
%plot(pCent,mN_Pw,'-b');hold on;plot(pCent,mN_Pint,'-r');
xlim([0 12000]);ylim([0 0.3]);
xlabel('P^{C} (atm)'); ylabel('{\it{f}}_{observed}');
yticks([0:0.1:0.2]);xticks([0:4000:12000]);
set(gca,'TickLength',[0.01 0.01],'LineWidth',2,'box','off','Layer','top');

%eval(['print -dtiff Fig1bc.tiff'])
%print(gcf,'-depsc2', '-painters', 'Fig1bc.eps');

%% Global otoliths  (Fig 5)
%For map
swloadpdf=1; %load existing data, or set to value other than one to regenerate (slow)
%Probability object vectors
dxw=0.2;bed=[(-25-dxw/2):dxw:(10+dxw/2)];bcent=[-25:dxw:10]; 
%change these only if similarly changed and reprocessed in the Global 
%calculation script

%Generate weighted average across depth
v=WOAsub.v3d;id=isnan(d13oto_3d_AEave);
if swloadpdf~=1
for i=1:size(WOAsub.temp,1) %longitude
    i
    for j=1:size(WOAsub.temp,2) %latitude
        clearvars binW binWint wtv ctZ;
        for k=1:size(WOAsub.temp,3) %depth    
            binW(k,:)=squeeze(wtPDFall(i,j,k,:));%raw d13 distribution weights
        end
        ctZ=nansum(binW,2).*(0.25.*0.2.*1)./(sum(Trt.tOMdist(:)*0.25)*sum(Trt.AEdist(:)*0.2.*1)); %weight by fraction of species present
        binWint=round(nansum(binW.*ctZ,1).*1e200,0); 
        %ensure integers along d13 bin vector
        %weighted by proportion of possible trait combinations contributing
        %to a particular WOA cell 
        if nansum(binWint)==0
            d13Z(i,j,1)=nan; 
        else
            pDw=fitdist(bcent','Kernel','Width',dxw,'Frequency',binWint'); 
            %pdf of combined d13 distribution weights in depth bin
            d13Z(i,j,:)=nanmean(random(pDw,1,5000)); 
            %generate distribution of d13oto values for each x and y, and
            %average across all depth cells at location
        end
    end
end
    save F5d13pdfZ_ct.mat d13Z -v7.3
else 
    load F5d13pdfZ_ct.mat
end

xx=double(WOAsub.x);xx=cat(1,xx(181:end),xx(1:180));
aa=d13Z';aax=cat(2,aa(:,181:end),aa(:,1:180));

f2=figure(2); 
 set(f2,'Renderer', 'painters', 'Position', [10 10 500 1100]);
 set(gcf,'color','w');
set(0,'DefaultAxesFontSize', 12);
set(0,'defaultLineLineWidth', 2);
%h=subplot(311);
h=subplot('Position',[0.08 0.68 0.9 0.32]);
axesm('MapProjection','hatano','Frame','on');
setm(h,'FLineWidth',1.5,'Grid','off');
tightmap;

surfm((double(WOAsub.y)-0.5),(xx-0.5),aax);hold on;

load coastlines;
geoshow(coastlat,coastlon,"DisplayType","polygon","FaceColor",'white');
cc=parula(12);colormap(gca,cc);cbh1=colorbar;
caxis([-6 0]);
%xlabel('Longitude (^{o}E)');ylabel('Latitude (^{o}N)');
ylabel(cbh1,'Mean \delta^{13}C_{oto}');
axis off;

clearvars xx aa aax d13Z i j;

%% Generate variation across traits and diet
%For state space diagram
swloadpdf=1; %if 1, load already calculated fields, else recalculate
t3d=nanmean(WOAsub.temp,4);
[~,~,bT]=histcounts(t3d(:),Goc.Tedge);
if swloadpdf~=1
for j=1:length(Goc.Tcent)%1:length(Goc.Tcent) %bx=j
    tic;j
    k=find(bT==j); %operate in one temperature bin
    length(k)
    if length(k)>0
        clearvars binW ctT binWint;
        for pix=1:length(k)%for each index location calculate distribution 
            [i1,i2,i3]=ind2sub(size(id),k(pix));
            binW(pix,:)=squeeze(wtPDFall(i1,i2,i3,:)); %d13 distribution weights
        end
        ctT=nansum(binW,2).*(0.25.*0.2.*1)./(sum(Trt.tOMdist(:)*0.25)*sum(Trt.AEdist(:)*0.2.*1));
        binWint=round(nansum(binW.*ctT,1).*1e200,0); 
        %ensure integers along d13 bin vector, within a single T bin, 
        %weighted by proportion of possible trait combinations contributing
        %to a particular WOA cell   
                    
        if nansum(binWint)==0
            d13T(j,1:5000)=nan; 
            K(j)=nan;
        else
            K(j)=length(k);
        pDw=fitdist(bcent','Kernel','Width',dxw,'Frequency',binWint'); 
        %pdf of combined volume-weighted d13 distribution weights in temperature bin
        d13T(j,:)=random(pDw,1,5000); %generate distribution of d13oto values for each temperature bin
        end
    else d13T(j,1:5000)=nan;binV(j)=nan;K(j)=nan; %no d13 values at all
    end
    %Create associated temperature vector
    Thb(j,:)=ones(1,5000).*Goc.Tcent(j);
    toc;
end
save F5d13pdf_ct.mat Thb d13T binV K -v7.3;
else
load F5d13pdf_ct.mat %F5d13pdf_ct.mat %weighted by fraction of trait combinations present
end

%Now have d13oto estimates, within each temperature bin, constant count in
%each temperature bin. Bin in 2d state space.
dxw2=dxw;bed2=[(-25-dxw2/2):dxw2:(10+dxw2/2)];bcent2=[-25:dxw2:10];
[N,~,~,bx,by]=histcounts2(Thb(:),d13T(:),Goc.Tedge,bed2); iN=(N==0);N(iN)=nan;
wtN=nan(size(N));
for j=1:length(Goc.Tcent) %bx=j
    k=find(bx==j);
     wtN(j,:)=N(j,:).*K(j)./nansum(K); %weight by relative number of solutions in each temperature bin
     %so temperatures with fewer solutions have a decreased relative
     %intensity.
end
wtN(iN)=nan;wtN(wtN==0)=nan;wtN=wtN./nansum(wtN,2); %sums to 1

h2=subplot('Position',[0.08 0.34 0.9 0.33]);
imagesc(Goc.Tcent,bcent2,(wtN'.*100),'AlphaData',~isnan((wtN'.*100))); %linear scale
set(h2,'YDir','normal');
colorbar;caxis([0 8]);shading flat;
xlabel('Longitude (^{o}E)');ylabel('Latitude (^{o}E)');
shading flat; hcb=colorbar;hcb.Ticks=[0:2:20];ylim([-8 4]);xlim([-2 30]); 
%ocean goes to -3C, but vanishingly small volume creates colorbar discontinuity
xlabel('Temperature (^{o}C)'); ylabel('\delta^{13}C_{oto}'); ylabel(hcb,'% species');
hold on;cc=parula(8);colormap(gca,cc);

clearvars d13T Thb N bx by;
%% Add observational data
%Otolith data overlay

for i=1:length(C_IS)
    if i==1 clr=[0.5 0.5 1];      %Adult Patagonian toothfish, wild
    elseif i==2 clr=[1 0 0];      %larval to juvenile Australian salmon, reared
    elseif i==3 clr=[1 0 1];      %larval red drum, reared
    elseif i==4 clr=[0 0 1];      %Juvenile and adult Pacific cod, wild
    elseif i==5 clr=[0 0 0];      %juvenile Atlantic cod, reared
    elseif i==6 clr=[0 0 0];      %larval to adult Atlantic cod, wild
    elseif i==7 clr=[0 0 0];      %adult Atlantic cod, reared
    elseif i==8 clr=[0 0 0];      %adult Atlantic cod, wild
    elseif i==9 clr=[0.8 0.8 1];  %adult orange roughy, wild
    elseif i==10 clr=[0.75 0.75 0.1]; %adult red emperor snapper, wild     
    elseif i==11 clr=[0.95 0.95 0.1]; %adult white-splotched grouper, wild
    elseif i==12 clr=[0.49 0.18 0.56]; %juvenile conger eel, wild
    elseif i==13 clr=[1 0.5 0.5]; %Pagrus auratus, wild
    elseif i==14 clr=[1 0.8 0.8]; %Pomatomus saltatrix, wild
    elseif i==15 clr=[0.7 0.7 0.7]; %Makaira nigricans, wild
    elseif i==16 clr=[1 1 1];     %Kajikia albidus, wild        
    end

lw=1.5; %line width for plot
ms=8; %marker size for plot

    %concatenate cod datasets
    if i==5
        tmp=cat(1,C_IS{5,1},C_IS{7,1});
    elseif i==6
        tmp=cat(1,C_IS{6,1},C_IS{8,1});
    elseif i==7 | i==8
        tmp=nan;
    else
        tmp=C_IS{i,1};
    end
    
    if all(isnan(tmp)) %do nothing
    else %plot
    t=tmp(:,2);m=tmp(:,3);o=tmp(:,1);
    %Reared or wild symbol
        if i==2|i==3|i==5; mkr='s'; %reared data 
        else mkr='o'; %wild data
        end
    %Fixed temperatures?
        if i==1|i==2|i==5|i==6
        [C,~,iC]=unique(t);
        for j=1:length(C)
            idx=iC==j;
            m=nanmean(o(idx));ll=nanmin(o(idx));ul=nanmax(o(idx));
        plot([C(j) C(j)],[ll,ul],'-','LineWidth',lw,'Color','k');hold on;
        plot(C(j),m,mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','k','MarkerSize',ms); 
        end
        end
    %Narrow temperature range, not already binned
        if i==4|i==9
            m=nanmean(o);ll=nanmin(o);ul=nanmax(o);
            tm=nanmean(t);tll=nanmin(t);tul=nanmax(t);
        plot([tm tm],[ll,ul],'-','LineWidth',lw,'Color','k');hold on;
        plot([tll tul],[m,m],'-','LineWidth',lw,'Color','k');
        plot(tm,m,mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','k','MarkerSize',ms); 
        end
    %Bespoke depending on peculiarities of dataset
        if i==3|i==12 %plot range for each bin
        for j=1:length(t)
        plot([t(j) t(j)],[tmp(j,6),tmp(j,7)],'-','LineWidth',lw,'Color','k');hold on;
        plot(t(j),o(j),mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','k','MarkerSize',ms);  
        end
        end
        if i==10|i==11 %group by T and then plot range for combined bins
        [C,~,iC]=unique(t);
        for j=1:length(C)
            idx=iC==j;
            m=nanmean(o(idx));ll=nanmin(nanmin(tmp(idx,6)));ul=nanmax(nanmax(tmp(idx,7)));
        plot([C(j) C(j)],[ll,ul],'-','LineWidth',lw,'Color','k');hold on;
        plot(C(j),m,mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','k','MarkerSize',ms); 
        end
        end
        if i==13|i==14 %group by T and then plot range for combined bins
            %but these datasets also have estuarine values that must be
            %removed
        flg=C_IS{i,5};tmp(flg==1,:)=[];t=tmp(:,2);m=tmp(:,3);o=tmp(:,1);
        [C,~,iC]=unique(t);
        for j=1:length(C)
            idx=iC==j;
            m=nanmean(o(idx));ll=nanmin(nanmin(tmp(idx,6)));ul=nanmax(nanmax(tmp(idx,7)));
        plot([C(j) C(j)],[ll,ul],'-','LineWidth',lw,'Color','k');hold on;
        plot(C(j),m,mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','k','MarkerSize',ms); 
        end
        end
        if i==15|i==16 %Marlin, only large temperature range known
            %Use commented out below if treating regions separately, but...
%         for j=1:length(t)
%             m=o(j);ll=(nanmin(tmp(j,6)));ul=(nanmax(tmp(j,7)));
%             tm=nanmean([tmp(j,10),tmp(j,11)]);tll=tmp(j,10);tul=tmp(j,11);
%         plot([tm tm],[ll,ul],'-','LineWidth',2,'Color','k');hold on;
%         plot([tll tul],[m,m],'-','LineWidth',2,'Color','k');
%         plot(tm,m,'o','LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','k','MarkerSize',6);  
%         end
%       ... here we bin all 3 regions together since no associated
%       temperature data and we have to estimate from regional
%       climatologies
            m=nanmean(o);ll=nanmin(o);ul=nanmax(o);
            tll=nanmin(tmp(:,10),[],'all');tul=nanmin(tmp(:,11),[],'all');
            tm=nanmean([tll,tul]);
        plot([tm tm],[ll,ul],'-','LineWidth',lw,'Color','k');hold on;
        plot([tll tul],[m,m],'-','LineWidth',lw,'Color','k');
        plot(tm,m,mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','k','MarkerSize',ms);
        end        
    end       
end %end plot loop

%The following is Kludgy, but...
%Repeat plot loop for mean data only, to avoid errorbars from one observation
%overlaying symbols from another
for i=1:length(C_IS)
    if i==1 clr=[0.5 0.5 1];          %Adult Patagonian toothfish, wild
    elseif i==2 clr=[1 0 0];      %larval to juvenile Australian salmon, reared
    elseif i==3 clr=[1 0 1];      %larval red drum, reared
    elseif i==4 clr=[0 0 1];      %Juvenile and adult Pacific cod, wild
    elseif i==5 clr=[0 0 0];%juvenile Atlantic cod, reared
    elseif i==6 clr=[0 0 0];%larval to adult Atlantic cod, wild
    elseif i==7 clr=[0 0 0];%adult Atlantic cod, reared
    elseif i==8 clr=[0 0 0];%adult Atlantic cod, wild
    elseif i==9 clr=[0.8 0.8 1];  %adult orange roughy, wild
    elseif i==10 clr=[0.75 0.75 0.1]; %adult red emperor snapper, wild     
    elseif i==11 clr=[0.95 0.95 0.1];%adult white-splotched grouper, wild
    elseif i==12 clr=[0.49 0.18 0.56]; %juvenile conger eel, wild
    elseif i==13 clr=[1 0.5 0.5]; %Pagrus auratus, wild
    elseif i==14 clr=[1 0.8 0.8]; %Pomatomus saltatrix, wild
    elseif i==15 clr=[0.7 0.7 0.7]; %Makaira nigricans, wild
    elseif i==16 clr=[1 1 1]; %Kajikia albidus, wild        
    end
    
    %concatenate cod datasets
    if i==5
        tmp=cat(1,C_IS{5,1},C_IS{7,1});
    elseif i==6
        tmp=cat(1,C_IS{6,1},C_IS{8,1});
    elseif i==7 | i==8
        tmp=nan;
    else
        tmp=C_IS{i,1};
    end
    
    if all(isnan(tmp)) %do nothing
    else %plot
    t=tmp(:,2);m=tmp(:,3);o=tmp(:,1);
    %Reared or wild symbol
        if i==2|i==3|i==5; mkr='s'; %reared data 
        else mkr='o'; %wild data
        end
    %Fixed temperatures?
        if i==1|i==2|i==5|i==6
        [C,~,iC]=unique(t);
        for j=1:length(C)
            idx=iC==j;
            m=nanmean(o(idx));ll=nanmin(o(idx));ul=nanmax(o(idx));
        plot(C(j),m,mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','k','MarkerSize',ms); 
        end
        end
    %Narrow temperature range, not already binned
        if i==4|i==9
            m=nanmean(o);ll=nanmin(o);ul=nanmax(o);
            tm=nanmean(t);tll=nanmin(t);tul=nanmax(t);
        plot(tm,m,mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','k','MarkerSize',ms); 
        end
    %Bespoke depending on peculiarities of dataset
        if i==3|i==12 %plot range for each bin
        for j=1:length(t)
        plot(t(j),o(j),mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','k','MarkerSize',ms);  
        end
        end
        if i==10|i==11 %group by T and then plot range for combined bins
        [C,~,iC]=unique(t);
        for j=1:length(C)
            idx=iC==j;
            m=nanmean(o(idx));ll=nanmin(nanmin(tmp(idx,6)));ul=nanmax(nanmax(tmp(idx,7)));
        plot(C(j),m,mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','k','MarkerSize',ms); 
        end
        end
        if i==13|i==14 %group by T and then plot range for combined bins
            %but these datasets also have estuarine values that must be
            %removed
        flg=C_IS{i,5};tmp(flg==1,:)=[];t=tmp(:,2);m=tmp(:,3);o=tmp(:,1);
        [C,~,iC]=unique(t);
        for j=1:length(C)
            idx=iC==j;
            m=nanmean(o(idx));ll=nanmin(nanmin(tmp(idx,6)));ul=nanmax(nanmax(tmp(idx,7)));
        plot(C(j),m,mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','k','MarkerSize',ms); 
        end
        end
        if i==15|i==16 %Marlin, only large temperature range known
            m=nanmean(o);ll=nanmin(o);ul=nanmax(o);
            tll=nanmin(tmp(:,10),[],'all');tul=nanmin(tmp(:,11),[],'all');
            tm=nanmean([tll,tul]);
        plot(tm,m,mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','k','MarkerSize',ms);
        end        
    end       
end %end plot loop

%% Chemistry and Eo lines
%diet for slopes
dint=[-16.85,-15.8,-15.8,-18.5,-18.5];
%slope intercepts
ax3=subplot('Position',[0.08 0.1 0.9 (0.557*0.33)]);
wtv=v./nansum(v(:));wtv=wtv(:);
%Arrhenius temperature term
Trange=[-3:1:30];kb=8.617e-5;
TT=(1/kb).*(1./(Trange+273.15)-1./(15+273.15));    
pO2ref=0.209;pCO2ref=400.*1e-6;
%Fish fluid ionic composition
pHf=(7.85-0.15);S_f=9;%Salinity equivalent of blood ionic strength for gas chemistry
[K0e K1e K2e]=carbeq(Trange',34);[K0f K1f K2f]=carbeq(Trange',S_f);
Khe=(O2sol(34,Trange)./(pO2ref).*1e-6)'; Khf=(O2sol(S_f,Trange)./(pO2ref).*1e-6)';%mol kg-1 atm-1 of O2 
[Dco2e,Dhco3e,Dco3e]=co2_diffusion(34,Trange');[Dco2f,Dhco3f,Dco3f]=co2_diffusion(S_f,Trange');
[~,~,~,~,Do2e,~,~]=gas_diffusion(34,Trange');[~,~,~,~,Do2f,~,~]=gas_diffusion(S_f,Trange');
qd=1*0.5;qs=K0e.*Dco2e./(Khe.*Do2e);
stn=0.011180;
%NOTE here that values of VhPhic are somewhat arbitrary for the excercise
%of comparing lines (though this and Eo could be fit for community using 
%supplemental figure). But here, just need to adjust intercepts for
%divergent lines and easiest way to do so is added offset to dmf.
for ci=1:5
    if ci==1
def=1.33;dmf=dint(1); %%%tweak here for plot, mean ~+1.33
Eo=0.35;%%%tweak here for plot, mean ~0.35%%%
VhPhic=0.286; %%%tweak here for plot, mean ~0.286%%%
clrl='k';
    elseif ci==2
def=1.33;dmf=dint(2);
Eo=0;%Eo=0 if chemistry only
VhPhic=0.286; 
clrl='r';
    elseif ci==3
def=1.33;dmf=dint(3); 
Eo=0;
VhPhic=0.286; 
clrl='r';
    elseif ci==4 %Just Eo and Vh*Phic for ocean <=0C
def=1.33;dmf=dint(4);
iTlow=t3d<=0; %subset to mean community traits in cells < or = 0C
Tlow=t3d;Tlow(~iTlow)=nan;
VhPhic_low=VhPhic_3d;VhPhic_low(~iTlow)=nan;
Eo_low=Eo_3d;Eo_low(~iTlow)=nan;
%generate weighted mean of thse traits
vpedge=[0.1:0.01:0.16];vpcent=vpedge(2:end)-mean(diff(vpedge))/2;
eoedge=[0.4:0.01:0.47];eocent=eoedge(2:end)-mean(diff(eoedge))/2;
[N2,~,~,bx2,by2]=histcounts2(Tlow(:),VhPhic_low(:),Goc.Tedge,vpedge); iN2=(N2==0);N2(iN2)=nan;
[N3,~,~,bx3,by3]=histcounts2(Tlow(:),Eo_low(:),Goc.Tedge,eoedge); iN3=(N3==0);N3(iN3)=nan;
wtN2=nan(size(N2));wtN3=nan(size(N3));
for i=1:length(vpcent) %by=i
    for j=1:length(Goc.Tcent) %bx=j
        k=find(by2==i&bx2==j);
        wtN2(j,i)=nansum(wtv(k));
    end
end
wtN2(iN2)=nan;wtN2(wtN2==0)=nan;
for i=1:length(eocent) %by=i
    for j=1:length(Goc.Tcent) %bx=j
        k=find(by3==i&bx3==j);
        wtN3(j,i)=nansum(wtv(k));
    end
end
wtN3(iN3)=nan;wtN3(wtN3==0)=nan;
VhPhic_lm=sum(nansum(wtN2,1)./nansum(wtN2,'all').*vpcent);
Eo_lm=sum(nansum(wtN3,1)./nansum(wtN3,'all').*eocent);
Eo=Eo_lm; %~0.42
VhPhic=VhPhic_lm; %~0.14
clrl='b';

    end


%Mass balance calculations
Ref=((def)./1000+1).*stn;Rmf=((dmf)./1000+1).*stn;
for i=1:length(Trange)%Ignoring small pH effect here and later
[~,ReCw,~,~,~,~]=C13speciation(Ref,8.1,K1e(i),K2e(i),Trange(i),'DIC'); %CO2aq,seawater
Rew(i,1)=ReCw;
end
PcritO2active=VhPhic.*exp(-Eo.*TT)'; %Vh(T)*PhiCrit
PcritC=PcritO2active.*(qd./qs);
Rint=nan(length(Trange),1);
d13_int=nan(length(Trange),1);
d13_sp=nan(length(Trange),1);
for tt=1:length(Trange)
Rint(tt,1)=(Rmf.*PcritC(tt)./pCO2ref+Rew(tt))./(PcritC(tt)./pCO2ref+1);%length(Ac)*length(Eo)*length(Rmf)
d13_int(tt,1)=(Rint(tt)./stn-1)*1e3;
[~,~,~,~,Rarag,~]=C13speciation(Rint(tt),pHf,K1f(tt),K2f(tt),Trange(tt),'CO2aq');
d13_sp(tt,1)=(Rarag./stn-1)*1e3; 
end
if ci==1
    plot(Trange,d13_sp,'-','Color',clrl,'LineWidth',1);hold on;
elseif ci==2
    plot(Trange(18:end),d13_sp(18:end),'-','Color',clrl,'LineWidth',1);hold on;
elseif ci==4
    plot(Trange(1:28),d13_sp(1:28),'-','Color',clrl,'LineWidth',1);hold on;
end
end

%Replot observation means only
    dscl=[-20.5:0.5:-16.5];
    cscl=repmat(linspace(0,0.87,length(dscl))',1,3);
for i=1:length(C_IS)    
    %concatenate cod datasets
    if i==5
        tmp=cat(1,C_IS{5,1},C_IS{7,1});
    elseif i==6
        tmp=cat(1,C_IS{6,1},C_IS{8,1});
    elseif i==7 | i==8
        tmp=nan;
    else
        tmp=C_IS{i,1};
    end
    
    if all(isnan(tmp)) %do nothing
    else %plot
    t=tmp(:,2);o=tmp(:,1);d=tmp(:,5);s=tmp(:,4);
    %Reared or wild symbol
        if i==2|i==3|i==5; mkr='s'; %reared data 
        else mkr='o'; %wild data
        end
    %Fixed temperatures?
        if i==1|i==2|i==5|i==6
        [C,~,iC]=unique(t);
        for j=1:length(C)
            idx=iC==j;
            mo=nanmean(o(idx));md=nanmean(d(idx));[~,icl]=min(abs(dscl-(md)));clr=cscl(icl,:);
            plot(C(j),mo,mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','none','MarkerSize',6);hold on;
        end
        end
    %Narrow temperature range, not already binned
        if i==4|i==9
            mo=nanmean(o);md=nanmean(d);tm=nanmean(t);[~,icl]=min(abs(dscl-(md)));clr=cscl(icl,:);
            plot(tm,mo,mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','none','MarkerSize',6);hold on;
        end
    %Bespoke depending on peculiarities of dataset
        if i==3|i==12 %plot range for each bin
        for j=1:length(t)
            [~,icl]=min(abs(dscl-(d(j))));clr=cscl(icl,:);
            plot(t(j),o(j),mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','none','MarkerSize',6);hold on;
        end
        end
        if i==10|i==11 %group by T and then plot range for combined bins
        [C,~,iC]=unique(t);
        for j=1:length(C)
            idx=iC==j;
            mo=nanmean(o(idx));md=nanmean(d(idx));[~,icl]=min(abs(dscl-(md)));clr=cscl(icl,:);
            plot(C(j),mo,mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','none','MarkerSize',6);hold on;
        end
        end
        if i==13|i==14 %group by T and then plot range for combined bins
            %but these datasets also have estuarine values that must be
            %removed
        flg=C_IS{i,5};tmp(flg==1,:)=[];t=tmp(:,2);m=tmp(:,3);o=tmp(:,1);
        [C,~,iC]=unique(t);
        for j=1:length(C)
            idx=iC==j;
            mo=nanmean(o(idx));md=nanmean(d(idx));[~,icl]=min(abs(dscl-(md)));clr=cscl(icl,:);
            plot(C(j),mo,mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','none','MarkerSize',6);hold on;      
        end
        end
        if i==15|i==16 %Marlin, only large temperature range known
            mo=nanmean(o);tll=nanmin(tmp(:,10),[],'all');tul=nanmin(tmp(:,11),[],'all');tm=nanmean([tll,tul]);
            md=nanmean(d);[~,icl]=min(abs(dscl-(md)));clr=cscl(icl,:);
            plot(tm,mo,mkr,'LineStyle','none','MarkerFaceColor',clr,'MarkerEdgeColor','none','MarkerSize',6);hold on;                   
        end        
    end 
end 
colormap(ax3,cscl);c=colorbar;c.Ticks=[dscl];

% labels
ylim([-6.7 0]);yticks([-6:2:0]);xlim([-2 30]);
xlabel('Temperature (^{o}C)'); ylabel('\delta^{13}C_{oto}');

set(findall(gcf,'-property','FontSize'),'FontSize',12)
% eval(['print -djpeg MI_otolith_Fig4_' figlab '.jpg'])

%% Maps (Fig 4 Suppl #1)
xx=double(WOAsub.x);xx=cat(1,xx(181:end),xx(1:180));

f4=figure(4);
 set(f4,'Renderer', 'painters', 'Position', [10 10 1000 (1.6*1000)]);
 set(gcf,'color','w');
set(0,'DefaultAxesFontSize', 12);
set(0,'defaultLineLineWidth', 0.5);

%correct for volume differences between hydrographic cells when plotting
%average concentrations
aa=nansum(dic_3d.*v./nansum(v,3),3)';aa(aa==0)=nan;aax=cat(2,aa(:,181:end),aa(:,1:180));
h41=subplot('Position',[0.01 0.55 0.42 0.45]);
axesm('MapProjection','hatano','Frame','on');
setm(h41,'FLineWidth',0.5,'Grid','off');
tightmap;surfm((double(WOAsub.y)-0.5),(xx-0.5),aax);hold on;
load coastlines;geoshow(coastlat,coastlon,"DisplayType","polygon","FaceColor",'white');
cc=parula(6);colormap(gca,cc);cbh41=colorbar;caxis([1900 2200]);
ylabel(cbh41,'DIC');axis off;clearvars aa aax;

aa=nansum(pco2_3d.*v./nansum(v,3),3)';aa(aa==0)=nan;aax=cat(2,aa(:,181:end),aa(:,1:180));
h42=subplot('Position',[0.51 0.55 0.42 0.45]);
axesm('MapProjection','hatano','Frame','on');
setm(h42,'FLineWidth',0.5,'Grid','off');
tightmap;surfm((double(WOAsub.y)-0.5),(xx-0.5),aax);hold on;
load coastlines;geoshow(coastlat,coastlon,"DisplayType","polygon","FaceColor",'white');
cc=parula(8);colormap(gca,cc);cbh42=colorbar;caxis([200 1000]);
ylabel(cbh42,'pCO_2');axis off;clearvars aa aax;

aa=nansum(d13dic_3d.*v./nansum(v,3),3)';aa(aa==0)=nan;aax=cat(2,aa(:,181:end),aa(:,1:180));
h43=subplot('Position',[0.01 0.05 0.42 0.45]);
axesm('MapProjection','hatano','Frame','on');
setm(h43,'FLineWidth',0.5,'Grid','off');
tightmap;surfm((double(WOAsub.y)-0.5),(xx-0.5),aax);hold on;
load coastlines;geoshow(coastlat,coastlon,"DisplayType","polygon","FaceColor",'white');
cc=parula(7);colormap(gca,cc);cbh43=colorbar;caxis([0.0 1.75]);set(cbh43,'XTick',[0:0.25:1.75]);
ylabel(cbh43,'\delta^{13}C_{DIC}');axis off;clearvars aa aax;

aa=nansum(d13co2_3d.*v./nansum(v,3),3)';aa(aa==0)=nan;aax=cat(2,aa(:,181:end),aa(:,1:180));
h44=subplot('Position',[0.51 0.05 0.42 0.45]);
axesm('MapProjection','hatano','Frame','on');
setm(h44,'FLineWidth',0.5,'Grid','off');
tightmap;surfm((double(WOAsub.y)-0.5),(xx-0.5),aax);hold on;
load coastlines;geoshow(coastlat,coastlon,"DisplayType","polygon","FaceColor",'white');
cc=parula(10);colormap(gca,cc);cbh44=colorbar;caxis([-12 -7]);
ylabel(cbh44,'\delta^{13}C_{CO2(aq)}');axis off;clearvars aa aax;

%eval(['print -djpeg MI_otolith_Fig4_Supp1.jpg'])

%% d13C components vs Latitude and Temp (Fig 4 Suppl #2)
Qs(Qs==0)=nan;
xl=[-70 70];

%169 to 170 degrees latitude
px=170;pk=[1];psym='-:'; %surface at 170 E

figure(5); clf reset; 
for kk=1:length(pk)
subplot(421)
plot(WOAsub.y,nanmean(VhPhic_3d(:,:,pk(kk)),1),['r' psym(kk)]); hold on; %shading flat;colorbar
plot(WOAsub.y,nanmean(PcritO2_3d(:,:,pk(kk)),1),['b' psym(kk)]); xlim(xl);
legend('V_h*\Phi_c','P_{crit}^{O2}(T)','Location','North');  ylabel('pO_2 (atm)')
subplot(423)
plot(WOAsub.y,1./nanmean(Qs(:,:,pk(kk)),1),['k' psym(kk)]); hold on; %shading flat;colorbar
legend('Q(T)'); ylabel('Ratio');xlim(xl);
subplot(425)
yyaxis left
plot(WOAsub.y,nanmean(PcritC_3d(:,:,pk(kk)),1).*1e6,['r' psym(kk)]); hold on; %shading flat;colorbar
ylim([0 3000]); ylabel('pCO_2 (\muatm)')
yyaxis right
plot(WOAsub.y,nanmean(pco2_3d(:,:,pk(kk)),1),['b' psym(kk)]); hold on; %shading flat;colorbar
legend('P_{met}^{CO2}','P_w^{CO2}');ylabel('pCO_2 (\muatm)')
ylim([0 3000]); xlim(xl); 
subplot(427)
plot(WOAsub.y,squeeze(nanmean(d13oto_3d_AEave(:,:,pk(kk)),1))','c'); xlim(xl); hold on;
plot(WOAsub.y,squeeze(nanmean(d13int_3d_AEave(:,:,pk(kk)),1))','r'); xlim(xl); hold on;
plot(WOAsub.y,squeeze(nanmean(d13co2_3d(:,:,pk(kk)),1))','b'); 
plot(WOAsub.y,squeeze(nanmean(d13dic_3d(:,:,pk(kk)),1))','b:'); 
plot(WOAsub.y,squeeze(nanmean(d13om_3d_DODave(:,:,pk(kk)),1))','-','Color',[0.92 0.6 0.1]); 
xlabel('Latitude');ylabel('\delta^{13}C');ylim([-22 5])
legend('Otolith','Internal','CO_2 (aq)','DIC','OM','Location','West')

subplot(422)
temp=nanmean(nanmean(WOAsub.temp(:,:,pk(kk),:),4),1); x2=[-3 30];treg=[-3:0.25:30];
[~,~,bTi]=histcounts(temp,treg);

clearvars aa bb cc;
aa=nanmean(VhPhic_3d(:,:,pk(kk)),1);bb=[aa',bTi'];
for xi=1:length(treg)
    ji=bb(:,2)==xi;
    cc(xi,1)=nanmean(bb(ji,1));
end
treg2=treg;treg2(isnan(cc))=[];cc(isnan(cc))=[];
plot(treg2,cc,['r' psym(kk)]); hold on; 

clearvars aa bb cc;
aa=nanmean(PcritO2_3d(:,:,pk(kk)),1);bb=[aa',bTi'];
for xi=1:length(treg)
    ji=bb(:,2)==xi;
    cc(xi,1)=nanmean(bb(ji,1));
end
treg2=treg;treg2(isnan(cc))=[];cc(isnan(cc))=[];
plot(treg2,cc,['b' psym(kk)]); xlim(xl);
legend('V_h*\Phi_c','P_{crit}^{O2}(T)','Location','West'); ylabel('pO_2 (atm)');xlim(x2);

subplot(424)
clearvars aa bb cc;
aa=1./nanmean(Qs(:,:,pk(kk)),1);bb=[aa',bTi'];
for xi=1:length(treg)
    ji=bb(:,2)==xi;
    cc(xi,1)=nanmean(bb(ji,1));
end
treg2=treg;treg2(isnan(cc))=[];cc(isnan(cc))=[];
plot(treg2,cc,['k' psym(kk)]); hold on;
legend('Q(T)','Location','NorthWest'); ylabel('Ratio');xlim(x2);

subplot(426)
clearvars aa bb cc;
aa=nanmean(PcritC_3d(:,:,pk(kk)),1)*1e6;bb=[aa',bTi'];
for xi=1:length(treg)
    ji=bb(:,2)==xi;
    cc(xi,1)=nanmean(bb(ji,1));
end
treg2=treg;treg2(isnan(cc))=[];cc(isnan(cc))=[];
yyaxis left
plot(treg2,cc,['r' psym(kk)]); hold on; 
clearvars aa bb cc;
aa=nanmean(pco2_3d(:,:,pk(kk)),1);bb=[aa',bTi'];
for xi=1:length(treg)
    ji=bb(:,2)==xi;
    cc(xi,1)=nanmean(bb(ji,1));
end
ylim([0 3000]); xlim(x2);
treg2=treg;treg2(isnan(cc))=[];cc(isnan(cc))=[];
yyaxis right
plot(treg2,cc,['b' psym(kk)]); hold on;
legend('P_{met}^{CO2}','P_w^{CO2}','Location','NorthWest');ylabel('pCO_2 (\muatm)');
ylim([0 3000]); xlim(x2);

subplot(428)
clearvars aa bb cc;
aa=nanmean(d13oto_3d_AEave(:,:,pk(kk)),1);bb=[aa',bTi'];
for xi=1:length(treg)
    ji=bb(:,2)==xi;
    cc(xi,1)=nanmean(bb(ji,1));
end
treg2=treg;treg2(isnan(cc))=[];cc(isnan(cc))=[];
plot(treg2,cc,['c' psym(kk)]); hold on;

clearvars aa bb cc;
aa=nanmean(d13int_3d_AEave(:,:,pk(kk)),1);bb=[aa',bTi'];
for xi=1:length(treg)
    ji=bb(:,2)==xi;
    cc(xi,1)=nanmean(bb(ji,1));
end
treg2=treg;treg2(isnan(cc))=[];cc(isnan(cc))=[];
plot(treg2,cc,['r' psym(kk)]); hold on;

clearvars aa bb cc;
aa=nanmean(d13co2_3d(:,:,pk(kk)),1);bb=[aa',bTi'];
for xi=1:length(treg)
    ji=bb(:,2)==xi;
    cc(xi,1)=nanmean(bb(ji,1));
end
treg2=treg;treg2(isnan(cc))=[];cc(isnan(cc))=[];
plot(treg2,cc,['b' psym(kk)]); hold on;

clearvars aa bb cc;
aa=nanmean(d13dic_3d(:,:,pk(kk)),1);bb=[aa',bTi'];
for xi=1:length(treg)
    ji=bb(:,2)==xi;
    cc(xi,1)=nanmean(bb(ji,1));
end
treg2=treg;treg2(isnan(cc))=[];cc(isnan(cc))=[];
plot(treg2,cc,['b:' ]); hold on;

clearvars aa bb cc;
aa=nanmean(d13om_3d_DODave(px,:,pk(kk)),1);bb=[aa',bTi'];
for xi=1:length(treg)
    ji=bb(:,2)==xi;
    cc(xi,1)=nanmean(bb(ji,1));
end
treg2=treg;treg2(isnan(cc))=[];cc(isnan(cc))=[];
plot(treg2,cc,['y' psym(kk)],'Color',[0.92 0.6 0.1]); hold on;
xlabel('Temperature (^{o}C)');ylabel('\delta^{13}C');ylim([-22 5]);xlim(x2)
end

%eval(['print -djpeg MI_otolith_Fig4_Supp2.jpg'])

return

