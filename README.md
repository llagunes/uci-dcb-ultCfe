# uci-dcb-ultCfe
Information for Ultrasensitivity and conformational free energy 

% =====================================================
General information: 
% =====================================================
Script files for MWC and non-allosteric model code are written in MatLab. All scripts and functions have comments to show what they do and how calculations are performed. 


% =====================================================
How to run: 
% =====================================================

Will need to download/clone/pull all files. 

Script files were written so that parameters are edited once. Script files to run are: 
MWC system model: UltrasensitivityVsCFE_Script_022020_MWC
Non-allosteric system model: UltrasensitivityVsCFE_Script_022020_INDEP

All other files are functions used to generate all the figures for investigating the effect of conformational free energy on ultrasensitivity in multisite biological processs. 


Open script file. For example, open UltrasensitivityVsCFE_Script_022020_MWC.
Then edit the following parameters: 

plotTypes - 1,2, or 3
% 1- only dose response curves
% 2- heatmaps (increasing c) 
% 3- scatter plots 

maxRuns - whole number > 0
% number of simulations to run where each will create a vector c - activation parameter

LVal - positive real number; preferably greater than 30 (smoother run) 
% L - basal activation 

randType - 0,1,2
% Randomization of c's: 0-do not randomize, 1-randomize, 2- exp slip

minDegC + maxDegC
% min and max magnitude of c possible. Note that 0<c<1 when choosing these degrees. 

minDegAlpha + maxDegAlpha
% min and max magnitude of alpha possible. Note that alpha>1 when choosing these degrees. Code assumes alpha_i=1 for all i.

allNVals - natural number greater than 1
% All the values of n to consider

Tcost - positive real number
% Maintenance cost per site


% =====================================================
NOTES: 
% =====================================================

- Depending on what type of plot is defined, script files will output 1+ figures. 
- When first downloading/cloning/etc, code will run for a small sample set. 
- Code is designed to work with MatLab! not C,C++, Python, R, etc. 

