function [ allKs,allvValsMat,allDoseResCurve,allHills  ] = INDEP_DoseResponsePlots( allNVals, maxRuns,dVal, rangV, rangKs,randType )
% This function plots all the scatterplots of H vs v
% The output is a list of structures for each value on containing the dose response,
% H, and values of v.

% Inputs:
% allNVals - all values on n (number of sites)
% maxRuns - number of runs/simulations for c/alpha
% dVal - d basal activation
% rangV - magnitudes of v (aka 10^rangV)
% rangK - maginutes of k (aka 10^rangK
% Htype - method of calculating ultrasensitivity
% randType - randomizing c or not
% plotTypes - which plots will be outputed

% Outline:
% 1. Find the dose response curve for each run/v/k
    % A. Set up values of v and k
    % B. Generate vector u (enzyme concentration)
    % C. Find the hill number associated to that curve
% 2. Plot each dose response curve on the same figure


%% Function Code
% Define paramters/variables
amtsNs=size(allNVals,2); % total number of sites
allUs=[0:0.05:1, 1:0.1:10, 10:2:60]; siUs=size(allUs,2);
% Energy formular DeltaG = -r*T*log(prod(c))
for nVal = 1:amtsNs
    n=allNVals(nVal); % number of sites
    % ====================================================
    % 1. Find the dose response curve for each run/c/alpha
    % A. Set up values of c and alpha
    minDegC= rangV(1); maxDegC= rangV(2); % min and max magnitude of c
    minDegK= rangKs(1); maxDegK = rangKs(2); % min and max magnitude of alpha
    rVks = minDegK + (maxDegK-minDegK).*rand(maxRuns,n); % randomly chose the degree's
    vectKs=10.^rVks; % alpha - modification rate parameter
    % generate c vector
    if randType ==1 % randomly chose c
        rVc = minDegC + (maxDegC-minDegC).*rand(maxRuns,n); % randomly chose the degree's
        vVal=10.^rVc; % log-uniform for values of ci
    elseif randType==0 % create grid of c=cBar values that span all magnitudes from input
        rVc=linspace(minDegC,maxDegC,maxRuns); cVal_pre=10.^rVc;vVal=repmat(cVal_pre',1,n);
    else % split total cfe evenly among n sites
        T=298;  % temperature in Kelvin
        r= 1.98/1000;  %  Boltzmann constant r in units J mol^(-1) Kelvin^(-1)
        maxDeltG=-r*T*log((10^minDegC)^n); minDeltG=-r*T*log((10^maxDegC)^n);
        % array of all delta g's
        allDeltGs=linspace(minDeltG,maxDeltG,maxRuns);cVal_pre=exp(-allDeltGs./(n*r*T));
        vVal=repmat(cVal_pre',1,n);
    end
    drSt=zeros(maxRuns,siUs); % each row is a dr curve
    allHs=zeros(maxRuns,1); % store hill number for each run/c/alpha
    for nRun=1:maxRuns
        vectVs=vVal(nRun,:);vectks_i=vectKs(nRun,:);fInf=prod(vectVs)/(dVal+prod(vectVs));
        for uInd=1:siUs
            u=allUs(uInd); drSt_pre=doseResponse_Indep( dVal,vectVs,vectks_i,n,u)./fInf;
            drSt(nRun,uInd)= drSt_pre;% normalized-store the dose response values
        end
        % C. Find the hill number associated to that curve
        hillNo1= calcHillFunc_Indep(vectVs,dVal,vectks_i,n);
        % Store H. Each column is H for c
        allHs(nRun)=hillNo1;
    end
    % Items to export for each value of n
    allKs.(sprintf('n_%d', n))=vectKs; % vector of alphas
    allvValsMat.(sprintf('n_%d', n))=vVal; % vector of c's
    allDoseResCurve.(sprintf('n_%d', n)) = drSt; % Dose response for each vector (matrix)
    allHills.(sprintf('n_%d', n)) = allHs;% Hill number associated to each dose response   
end
%% Plot 
% ============================================
% 2. Plot each dose response curve on the same figure

% figure(plotTypes); clf; 
for nVal = 1:amtsNs
    n=allNVals(nVal);
    figure(10+n); clf;
    HsPlot=allHills.(sprintf('n_%d', n)); % extract hill number
    drCurveN=allDoseResCurve.(sprintf('n_%d', n)); % extract dr values
    for nRun=1:maxRuns        
        plot(allUs,drCurveN(nRun,:),'LineWidth',2); hold on; grid on;
        Hval=HsPlot(nRun);legendInfo{nRun} = ['H=' num2str(Hval,3)]; 
    end
    ylim([0,1.2]);
    xlabel('u: Enzyme concentration'); ylabel('f(u,c,\alpha): Dose response');
    hT=title([ 'Indep DR n= ' num2str(n) 'with d= ' num2str(dVal) '  k_i=1' ]);
    legend(legendInfo);set(gca,'fontsize',25,'fontWeight','bold');   
end








end


