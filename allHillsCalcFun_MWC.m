function [ Cs_ord_c , Hplot_cs, cVal ,allHs, vectAlphas] = allHillsCalcFun_MWC( n , maxRuns, LVal, rangC, rangAlpha, Htype,randType)
% This function returns all the Hill numbers for given a value of n along with
% the values of c and cBar.

% Outline:
% 1. Set up values of c and alpha
% 2. Calculate arithmetic mean cBar (and aBar if needed)
% 3. Calculate H
% 3 a. for each vector c for all the runs/pts in simulation
% 3 b. for vector c with c_i = cBar
% 4. Output all H and c,alpha, cBar values

% Htype = type of Hill number/ultrasensitivity calculation:
% Htype = 1 - Golbeter-Kohsland
% Htype = 2 - Hill function Fit
% Htype = 3 - max log slope
% Htype = 4 - Levitsky n_50

% randType = randomization of all c's
% randType= 0- Do not randomize, create grid along magnitudes of c
% randType= 1- Randomize all values of c log-uniformly across magnitudes
% randType= 2- Do not randomize, split total cfe evenly
% ====================================================
% 1. Set up values of c and alpha
minDegC= rangC(1); maxDegC= rangC(2); % min and max magnitude of c
minDegAlpha= rangAlpha(1); maxDegAlpha = rangAlpha(2); % min and max magnitude of alpha
rValpha = minDegAlpha + (maxDegAlpha-minDegAlpha).*rand(maxRuns,n); % randomly chose the degree's
vectAlphas=10.^rValpha; % alpha - modification rate parameter
% generate c vector
if randType ==1 % randomly chose c
    rVc = minDegC + (maxDegC-minDegC).*rand(maxRuns,n); % randomly chose the degree's
    cVal=10.^rVc; % log-uniform for values of ci
elseif randType==0 % create grid of c=cBar values that span all magnitudes from input
    rVc=linspace(minDegC,maxDegC,maxRuns); cVal_pre=10.^rVc;
    cVal=repmat(cVal_pre',1,n);
else % split total cfe evenly among n sites
    T=298;  % temperature in Kelvin
    r= 1.98/1000;  %  Boltzmann constant r in units J mol^(-1) Kelvin^(-1)
    maxDeltG=-r*T*log((10^minDegC)^n); minDeltG=-r*T*log((10^maxDegC)^n);
    % array of all delta g's
    allDeltGs=linspace(minDeltG,maxDeltG,maxRuns);cVal_pre=exp(-allDeltGs./(n*r*T));
    cVal=repmat(cVal_pre',1,n);
end
% ====================================================
% 2. Calculate arithmetic mean cBar (and aBar if needed)
% Matrix that contains each pair of c_i=cBar
cBar=mean(cVal,2);
% Store all the H
allHs=zeros(maxRuns,1);
% ====================================================
% 3. Calculate H
for runN=1:maxRuns % For each epsilon vector
    % Define the vector c,and alpha
    vectCs=cVal(runN,:); vectAlphas_i=vectAlphas(runN,:);
    % 3a. Calculate H_1 for each vector c for all the runs/pts in simulation
    % calculation depends of the Htype
    if Htype==1 % Goldbeter-Kohsland
        hillNo1= calcHillFunc(vectCs,LVal,n,vectAlphas_i);
    elseif Htype==2 % Fit to Hill function
        hillNo1= calcHillFunc_Fitted(vectCs,LVal,n,vectAlphas_i);
    elseif Htype==3 % max-log slope
        hillNo1= calcHillFunc_maxLogSlope(vectCs,LVal,n,vectAlphas_i);
    else % Levitsky n_H defition (Eq. 91)
        hillNo1= calcHillFunc_Levitzki(vectCs,LVal,n,vectAlphas_i);
    end
    % Store H. Each column is H for c
    allHs(runN)=hillNo1;
end
% end
% ====================================================
% 4. Output all H and c,alpha, cBar values
% For plotting purposes, need to sort vectors according to arithmetic mean
[Cs_ord, indCs] = sort(cBar, 'ascend');
H_ord_cs=allHs(indCs);

% Need to remove all points where H is undefined, aka H=0 or NaN
a=(H_ord_cs==0);
% aNaN=isnan(H_ord_cs)==1; bNaN=isnan(H_ord_cBar)==1;
Cs_ord_c=Cs_ord(a==0);Hplot_cs=H_ord_cs(a==0);


% Last Edit: 01/29/2020 LL
end

