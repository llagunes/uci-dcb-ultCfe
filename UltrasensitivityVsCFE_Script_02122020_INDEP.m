%% INDEPENDENT and ultrasensitivity 

% This file will generate figures in corresponding paper in which we have
% scatterplots of ultrasensitivity (H) vs the total conformational free
% energy (cfe) in the generalized non-allosteric system (independent)

%% NOTES:

% General outline:
% 1. Fix n and L and types of plots (Notes below)
% 2. Generate grid of values of c and alpha (if changing simultaneously)
% 3. Calculate arithmetic mean of c cBar (and aBar if needed) 
% 4. Calculate H_1 for each vector c for all the runs/pts in simulation
% 5. Calcualte H_2 for vector c with c_i = cBar
% 6. Plot H vs totalCFE

%% Code 
clear; clc;
% To store data
ModelType='INDEP'; HH=0212; % date - no dashes
allHTypes = {'H_{GK}', 'H_{Fit}', 'H_{MLS}' ,'H_{Lev}' };
Htype=1; % Type of H calculation: 1 - GoldbeterKoshland, 2- Fit to Hill, 3- max log slope, 4-Levitzki
% =======================================================
% 1. Fix n and L and types of plots
% Types of plot to ouput 
% 1- only dose response curves, 2- heatmaps (increasing c) , 3- scatter plots 
plotTypes=3; 
maxRuns = 100; % number of simulations to run where each will randomize c
dVal=1000; % L - basal activation 
% Parameters to set in conversion to energy units from the Gibbs Free
% Energy formular DeltaG = -r*T*log(prod(c))
T=298;  % temperature in Kelvin
r= 1.98/1000;  %  Boltzmann constant r in units J mol^(-1) Kelvin^(-1)
% =======================================================
% 2. a. Generate grid of values of c
% Here, each c_i will be log-uniformily drawn 
randType=0;% Randomization of c's: 0-do not randomize, 1-randomize, 2- exp slip
minDegV=1;maxDegV=4;% min and max magnitude of c possible
rangV=[minDegV, maxDegV];
% NOTE: max degree I could hit was -0.01 to get close to c=1.
% The case where both c1=c2=1 is not possible since that would mean no site
% contributes to substrate activation. So max degree is in [-4,0) although
% theoretically it can be (-oo,0).
% =======================================================
% 2. b. Generate grid of values of alpha, alpha - modification rate parameter
% Here, each alpha_i will be the same 
minDegK = 0;maxDegK=0;rangKs=[minDegK, maxDegK];
% For each value of n, calculate H
allNVals=[2,3,4,6,8];

if plotTypes==1 % Dose response curves only 
    [ allKs,allvValsMat,allDoseResCurve,allHills  ] = INDEP_DoseResponsePlots( allNVals, maxRuns,dVal, rangV, rangKs,randType );

elseif plotTypes == 2 % Heatmaps
    [ allvValsMat, allHsLog ] = INDEP_Heatmaps( allNVals, maxRuns,dVal, rangV, rangKs );

else  % Scatter plots
    [ allcValsMat,allvBARs,allCVs, allHsInd,allHsSame,allKsMat] = INDEP_ScatterPlots( allNVals, maxRuns,dVal, rangV, rangKs,plotTypes );
end 

%% END OF SCRIPT
% Last Edit: 02/12/2020 LL