function [ hillNo, maxLogSlope ] = calcHillFunc_maxLogSlope( vectCs,LVal,n,vectAlphas_i )
% This function calculates the Hill coefficient given the parameters above.
% Here, the Hill number is a measurement of ultrasensitivity found as the
% maximum log slope in a dose response curve. 

% OUTLINE: 
% 1. Find the Dose Response Function and ratio  output/input
% 2. Find log of DR ratios
% 3. Calculate the slope at each point 
% 4. Find the max slope

% ==========================================
% 1. Find the Dose Response Function and ratio  output/input
% f^Infinity = the max f value
fInf=1/(1+LVal*prod(vectCs));
% All values of u = inputs (aka enzyme concentrations) 
allUs = [0.01:.1:1, 2:5:20, 30:10:100, 200:100:1000, 2000:1000:10000]; 
maxUs=length(allUs); % Total number of points
% f evaluated at each u 
fDoseResponse=zeros(1,maxUs);
for indU = 1:maxUs
    u=allUs(indU);
    fDoseResponse(indU)=doseResponse_MWC( LVal,n,vectCs,vectAlphas_i,u);
end 
% ==========================================
% 2. Find log of DR ratios
% All values of f normalized to f^Infinity
f_norm = fDoseResponse./fInf;
ratioF=f_norm./allUs;% ratio of output/input
logRatioF=log(ratioF); 
% ==========================================
% 3. Calculate the slope at each point
SlopeBetweenPoints = diff(logRatioF)./diff(allUs);
% ==========================================
% 4. What is the max slope
maxLogSlope = max(SlopeBetweenPoints);
hillNo=exp(maxLogSlope);

% Last edit: 12/11/19
end

