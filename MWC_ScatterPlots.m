function [ allcValsMat,allcBARs, allHs, allTCFEs,allHverif,allAlphasMat,allKnees] = MWC_ScatterPlots( allNVals, maxRuns,LVal, rangC, rangAlpha, Htype,randType,plotTypes,Tcost,mtdi )
% This function plots all the scatterplots of H vs total cfe + a MaintCost
% The output is a list of structures for each value on containing totalCFE,
% H, and values of c.


% Inputs: 
% allNVals - all values on n (number of sites) 
% maxRuns - number of runs/simulations for c/alpha
% LVal - L basal activation 
% rangC - magnitudes of c (aka 10^rangC)
% rangAlpha - maginutes of alpha (aka 10^rangAlpha  
% Htype - method of calculating ultrasensitivity 
% randType - randomizing c or not 
% plotTypes - which plots will be outputed 
% Tcost - maintenance cost per site 
% mtdi - method of calculating the knee of a curve


% Outline: 
% 1: for each value of n, find H vs total cfe 
% 2: Plot figure

%% Function Code
% Define paramters/variables
amtsNs=size(allNVals,2); % total number of sites
% Energy formular DeltaG = -r*T*log(prod(c))
T=298;  % temperature in Kelvin
r= 1.98/1000;  %  Boltzmann constant r in units J mol^(-1) Kelvin^(-1)
for nVal = 1:amtsNs
    n=allNVals(nVal); % number of sites
    % ============================================
    % 1: for each value of n, find H vs total cfe    
    [Cs_ord_c_2 , Hplot_cs_2, cVals_2 ,allHs_2, vectAlphas_2]=allHillsCalcFun_MWC( n , maxRuns, LVal, rangC, rangAlpha, Htype,randType);
    % Also track the total conformational free energy
    totalCFE_n2_Hgk_pre=-r*T*log(prod(cVals_2,2));H_GK_all_n2=allHs_2(allHs_2>0);
    totalCFE_n2_Hgk=totalCFE_n2_Hgk_pre(allHs_2>0);% Remove undefined values of H
    % Store the key information in structures for each value of n
    allcValsMat.(sprintf('n_%d', n))=cVals_2;
    allcBARs.(sprintf('n_%d', n)) = Cs_ord_c_2;
    allHs.(sprintf('n_%d', n))=allHs_2;
    allTCFEs.(sprintf('n_%d', n))=totalCFE_n2_Hgk;
    allHverif.(sprintf('n_%d', n))=H_GK_all_n2;  
    allAlphasMat.(sprintf('n_%d', n))=vectAlphas_2;
    % ================================================================
    % Find where the curve begins to level off (aka the knee)
    % mtdi=3; % Sapotta 2011 method for finding the knee
    cutOffval=0.5;% Top % of points
    [ tableKnees ] = knees_MWC( totalCFE_n2_Hgk,H_GK_all_n2,n,mtdi,cutOffval );
    % tableKnees(1) = value of n,  tableKnees(2)= xvalue of where knee is, tableKnees(3)= Yvalue of knee,
    allKnees.(sprintf('n_%d', n))=tableKnees;
   
end 
%% Plot 
% ============================================
% 2: Plot figure

% 6. Make scatter plot with maint cost included 
% Plot H vs cBar + MAINT COST (Mc)
% Tcost=2; % Maint cost per site
colorsPlot={'k*','c-','g-','r-','k-'};
figure(plotTypes); clf; 
for nVal = 1:amtsNs
    n=allNVals(nVal);
    cVals=allcBARs.(sprintf('n_%d', n)); HsPlot=allHverif.(sprintf('n_%d', n));
    totalCFE_cBar_mc=-r*T*log(cVals.^n) + Tcost*n;
    % Need to extend the plots out further in their saturation level
    totalCFE_extd=[max(totalCFE_cBar_mc),1.15*max(totalCFE_cBar_mc)];
    if nVal==1
        plot(totalCFE_cBar_mc, HsPlot, '-', 'Color', [1, 0.6, 0.0],'LineWidth',2); hold on; grid on;
        % Plot the extension 
        plot(totalCFE_extd,[max(HsPlot),max(HsPlot)],'-', 'Color', [1, 0.6, 0.0],'LineWidth',2)
    else 
        plot(totalCFE_cBar_mc, HsPlot, colorsPlot{nVal},'LineWidth',2); hold on; grid on; 
        plot(totalCFE_extd,[max(HsPlot),max(HsPlot)], colorsPlot{nVal},'LineWidth',2)
    end 
    
end

%Ste5 + Mc 
plot(11+Tcost*8,4,'k*', 'LineWidth',3); hold on; % Ste5 data point

ylim([0.01,5.1]); %xlim([10^-4,10^-1])
xlabel('Total Conformational Free Energy'); ylabel('H');
hT=title([ 'H with L= ' num2str(LVal) '  \alpha_i=1 + Mc=' num2str(Tcost,1)]);
% set(hT,'interpreter','Latex','FontSize',25)
set(gca,'fontsize',25,'fontWeight','bold');

% Labels on Figure:
% Label H type
% x1 = 1;y1 = 0.5;txt1 = allHTypes(Htype);text(x1,y1,txt1,'fontsize',25,'Color','k','FontWeight','bold')
% n=2
x1 = 2;y1 = 1.5;% Mc=2
% x1 = 2;y1 = 1.5;% Mc=4
% x1 = 5;y1 = 1.5;% Mc=8
txt1 = 'n=2';text(x1,y1,txt1,'fontsize',25,'Color', [1, 0.6, 0.0],'FontWeight','bold')
% n=3
x1 = 6.8;y1 = 2.2;%Mc=2
% x1 = 10;y1 = 2.2;%Mc=4
% x1 = 15;y1 = 2.2;%Mc=8
txt1 = 'n=3';text(x1,y1,txt1,'fontsize',25,'Color','c','FontWeight','bold')
% n=4
x1 = 10;y1 = 2.9;% Mc=2
% x1 = 14;y1 = 2.9;% Mc=4
% x1 = 25;y1 = 2.9;% Mc=8
txt1 = 'n=4';text(x1,y1,txt1,'fontsize',25,'Color','g','FontWeight','bold')
% n=6
x1 = 17;y1 = 3.8;%Mc=2
% x1 = 28;y1 = 3.8;%Mc=4
% x1 = 45;y1 = 3.8;%Mc=8
txt1 = 'n=6';text(x1,y1,txt1,'fontsize',25,'Color','r','FontWeight','bold')
% n=8
x1 = 28;y1 = 4.53;%MC=2
% x1 = 50;y1 = 4.7;%MC=4
% x1 = 70;y1 = 4.7;%MC=8
txt1 = 'n=8';text(x1,y1,txt1,'fontsize',25,'Color','k','FontWeight','bold')
% STE5
x1 = 22;y1 = 4.2;%Mc=2
% x1 = 35;y1 = 4.2;%Mc=4
% x1 = 60;y1 = 4.2;%Mc=8
txt1 = 'Ste5';text(x1,y1,txt1,'fontsize',20,'Color','k','FontWeight','bold')





end

