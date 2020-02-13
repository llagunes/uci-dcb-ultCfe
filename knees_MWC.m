function [ tableKnees ] = knees_MWC( xVals, yVals,n,mtdi,cutOffval )
% This function will calculate the knee of the given curve from the input
% and output based on mtdi: 1-slopes, 2-curvature, 3-Sapotta. Output table
% format: Col1-n, col2-xValueKnee, col3-yValueKnee

% Knee calculations are: 
% 1- comparing slopes
% 2- maximum curvature
% 3- Sapotta paper

% OUTLINE: 
% 1. Find the knee (call function) when using only the top % of data set
% 2. Export table with knees

% ============================================================ 
% 1. Find the knee (call function) when using only the top % of data set
% with Sapotta & 80%
% cutOffval=0.8; 
top=1;
% mtdi=3;
[ x_p2_sap_t, y_knee2_80_sap_t] = calc_MaxCurv_Kneedle( xVals, yVals,cutOffval,mtdi,top);

% ============================================================ 
% 3. Export table with knees
% Order in table: 
% Col1 - value on n (need to have dimensions meet, hence the repetition)
% Col2 - x value of where knee occurs 
tableKnees = [n,x_p2_sap_t,y_knee2_80_sap_t];



% Last edit: 01/27/20 LL
end

