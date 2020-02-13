function [ allcValsMat,allvBARs,allCVs, allHsInd,allHsSame,allKsMat] = INDEP_ScatterPlots( allNVals, maxRuns,dVal, rangV, rangKs,plotTypes )
% This function plots all the scatterplots of H vs mean(v) 
% The output is a list of structures for each value on containing mean(v),
% H, cv(v) and values of v.


% Inputs: 
% allNVals - all values on n (number of sites) 
% maxRuns - number of runs/simulations for c/alpha
% dVal - d basal activation 
% rangV - magnitudes of c (aka 10^rangC)
% rangKs - maginutes of alpha (aka 10^rangAlpha  
% plotTypes - which plots will be outputed 

% Outline: 
% 1: for each value of n, find H vs total cfe 
% 2: Plot figure

%% Function Code
% Define paramters/variables
amtsNs=size(allNVals,2); % total number of sites

for nVal = 1:amtsNs
    n=allNVals(nVal); % number of sites
    % ============================================
    % 1: for each value of n, find H vs total cfe 
    [ actualEpsVec,meanEps,allCvs,hillsOrd,hillsOrdSame,vectAllKs ] = hillsVsMeanCV_INDEP( n, maxRuns, rangV,rangKs,dVal );
    % Store the key information in structures for each value of n
    allcValsMat.(sprintf('n_%d', n))=actualEpsVec;
    allvBARs.(sprintf('n_%d', n)) = meanEps;
    allCVs.(sprintf('n_%d', n)) = allCvs;
    allHsInd.(sprintf('n_%d', n))=hillsOrd;
    allHsSame.(sprintf('n_%d', n))=hillsOrdSame;
    allKsMat.(sprintf('n_%d', n))=vectAllKs;   
end 
%% Plot 
% ============================================
% 2: Plot figure

% 6. Make scatter plot with maint cost included 
% Plot H vs cBar + MAINT COST (Mc)
% Tcost=2; % Maint cost per site
colorsPlot={'k*','c.','g.','r.','k.'};
figure(plotTypes); clf; 
for nVal = 1:amtsNs
    n=allNVals(nVal);
    vVals=allvBARs.(sprintf('n_%d', n)); HsPlot=allHsInd.(sprintf('n_%d', n));
    HsPlot_Same=allHsSame.(sprintf('n_%d', n));
    if nVal==1
        semilogx(vVals, HsPlot, '.', 'Color', [1, 0.6, 0.0],'LineWidth',2); hold on; grid on;
        % Plot magenta line
        semilogx(vVals, HsPlot_Same, '-m', 'LineWidth',2); hold on; grid on;
    else
        semilogx(vVals, HsPlot_Same, '-m', 'LineWidth',2); hold on; grid on;
        semilogx(vVals, HsPlot, colorsPlot{nVal},'LineWidth',2); hold on; grid on;
    end
    
end

ylim([0.9,2])
xlabel('Mean(v)'); ylabel('H');
title(['H with d= ' num2str(dVal), ' and k_i=k_j '])
set(gca,'fontsize',25);
xlim([10^1,10^4])
% Add text to plot
% n=2
x1 = 10^3;y1 = 0.95;txt1 = 'n=2';text(x1,y1,txt1,'fontsize',25,'Color', [1, 0.6, 0.0],'FontWeight','bold')
% n=3
x1 = 10^3.6;y1 = 1.18;txt1 = 'n=3';text(x1,y1,txt1,'fontsize',25,'Color','c','FontWeight','bold')
% n=4
x1 = 10^3.6;y1 = 1.25;txt1 = 'n=4';text(x1,y1,txt1,'fontsize',25,'Color','g','FontWeight','bold')
% n=6
x1 = 10^3.6;y1 = 1.3;txt1 = 'n=6';text(x1,y1,txt1,'fontsize',25,'Color','r','FontWeight','bold')
% n=8
x1 = 10^3.6;y1 = 1.4;txt1 = 'n=8';text(x1,y1,txt1,'fontsize',25,'Color','k','FontWeight','bold')


end

