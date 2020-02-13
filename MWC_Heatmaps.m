function [ allcValsMat, allHsLog,allHsLin ] = MWC_Heatmaps( allNVals, maxRuns,LVal, rangC, rangAlpha )
% This function plots all the heatmaps of H vs c_1 and c_2 (by increasing
% only c_1 or c_2
% The output is a list of structures for each value on containing
% H and values of c.

% Inputs:
% allNVals - all values on n (number of sites)
% maxRuns - number of runs/simulations for c/alpha
% LVal - L basal activation
% rangC - magnitudes of c (aka 10^rangC)
% rangAlpha - maginutes of alpha (aka 10^rangAlpha
% Htype - method of calculating ultrasensitivity
% plotTypes - which plots will be outputed


% Outline:
% 1: Fix n, alpha, and L (and ci for i>=3)
% 2: Generate grid of values of c
% 3: Calculate H for each vector c for all the runs/pts needed
% 4: Plot H vs c1 and c2

%% Function Code
% ==========================================
% 1: Fix n, alpha, and L (and ci for i>=3)
amtsNs=size(allNVals,2); % total number of sites
eiVal=0.01; % c value for ci for n>2
minDegC= rangC(1); maxDegC= rangC(2); % min and max magnitude of c
% ==========================================
% 2: Generate grid of values of c
% NOTE: max degree I could hit was -0.01 to get close to c=1.
% The case where both c1=c2=1 is not possible since that would mean no site
% contributes to substrate activation. So max degree is in [-4,0) although
% theoretically it can be (-oo,0).
e1Vec=10.^(linspace(minDegC, maxDegC, maxRuns)); % all values of c along heatmap axis
e1VecLin= linspace(10^minDegC, 10^maxDegC,maxRuns); % in a linear scale
allcValsMat=[e1Vec',e1VecLin'];
% vector of alpha values
minDegAlpha= rangAlpha(1); maxDegAlpha = rangAlpha(2); % min and max magnitude of alpha
% for each value of n find H for each small increase on c_1/c_2 
for nVal=1:amtsNs 
    n=allNVals(nVal);   
    % Store all the H
    allHs=zeros(maxRuns,maxRuns);allHsLin_pre=zeros(maxRuns,maxRuns);
    rValpha = minDegAlpha + (maxDegAlpha-minDegAlpha).*rand(maxRuns,n); % randomly chose the degree's
    vectAlphas=10.^rValpha; % alpha - modification rate parameter
    % ==========================================
    % 3: Calculate H for each vector c for all the runs/pts needed
    for e2=1:maxRuns % For each c vector
        for e1=1:maxRuns % change c_2 first
            % For each c on the grid, make vector
            vectEPS=[e1Vec(e1),e1Vec(e2), eiVal*ones(1,n-2)];
            vectEPSLIN=[e1VecLin(e1),e1VecLin(e2), eiVal*ones(1,n-2)];
            % Calculate H for c_1 and c_2
            hillNo1= calcHillFunc(vectEPS,LVal,n,vectAlphas);
            [hillNoLin]= calcHillFunc(vectEPSLIN,LVal,n,vectAlphas);
            if isreal(hillNo1) == 0 % Remove imaginary values
                hillNo1=0;
            end
            % Store H. Each column is H for e1, row-e2.
            allHs(e2,e1)=hillNo1;allHsLin_pre(e2,e1)=hillNoLin;
        end
    end
    % Store H for this n (1 heatmap)
    allHsLog.(sprintf('n_%d', n))=allHs;allHsLin.(sprintf('n_%d', n))=allHsLin_pre;
end



%% Plot
% ===================================
% 4: Plot H vs c1 and c2

% ticks to label on heatmap
x1=1;x2=ceil((maxRuns)/2);x3=maxRuns;
%LOG
%All H=0 need to be NaN instead
for nVal=1:amtsNs 
    n=allNVals(nVal);
    matHs=flipud(allHsLog.(sprintf('n_%d', n)));
    mtHColor=matHs;mtHColor(mtHColor==0)= NaN; % Replace all 0's with NaN in H's matrix
    % numElmnts=sum(sum(mtHColor==0))/(maxRuns^2); %percent of elements that are NaN/0
    totNumElms=maxRuns^2;maxHVal=max(max(matHs));minHVal=min(min(mtHColor));
    
    figure(100+n); clf;
    im1=imagesc(mtHColor);
    colormap( [1 1 1; parula(totNumElms)] );colorbar;caxis( [minHVal/100 maxHVal] )
    yticklabels={'10^{-4}','10^{-2}','0.9'}; % 0.9 - 0.1
    yticks= [x1,x2,x3];set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
    xticklabels = yticklabels;
    xticks = [x1,x2,x3];set(gca, 'XTick', xticks, 'XTickLabel', (xticklabels(:)))
    % set(gca,'XScale','log'); set(gca,'YScale','log');
    ylabel('c_2'); xlabel('c_1')
    title(['H for n= ' num2str(n) ', L= ' num2str(LVal) ' \alpha_i=\alpha_j'] )
    set(gca,'FontSize',25,'fontWeight','bold')

    %LINEAR
    matHsLin=flipud(allHsLin.(sprintf('n_%d', n)));
    mtHColor=matHsLin;mtHColor(mtHColor==0)= NaN; % Replace all 0's with NaN in H's matrix
    maxHVal=max(max(matHsLin));minHVal=min(min(mtHColor));
    
    figure(200+n);clf;
    im1=imagesc(mtHColor );colormap( [1 1 1; parula(totNumElms)] )
    colorbar;caxis( [minHVal/100 maxHVal] )
    
    yticklabels={'10^{-4}','0.5','0.9'}; % 0.9 - 0.1
    yticks= [x1,x2,x3];set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
    xticklabels = yticklabels;
    xticks = [x1,x2,x3];set(gca, 'XTick', xticks, 'XTickLabel', (xticklabels(:)))
    % set(gca,'XScale','log'); set(gca,'YScale','log');
    ylabel('c_2'); xlabel('c_1')
    title(['H for n= ' num2str(n) ', L= ' num2str(LVal) ' \alpha_i=\alpha_j '])
    set(gca,'FontSize',25,'fontWeight','bold')
    
end 












end





