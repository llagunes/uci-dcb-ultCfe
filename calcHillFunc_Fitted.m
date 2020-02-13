function [ hillNo,ec50 ] = calcHillFunc_Fitted(vectCs,LVal,n,vectAlphas_i)
% This function will return the hill number by fitting to a non linear
% curve

% Outline: 
% 1. Generate the dose response curve with parameters c,alpha,L
% 2. Normalize curve
% 3. Fit to Hill Equation: x^n/(ec50^n+x^n), where n=Hill number

% =======================================
% 1. Generate the dose response curve with parameters c,alpha,L
% f^Infinity or the max f value
fInf=1/(1+LVal*prod(vectCs));
% All values of u to consider 
allUs = [0.01:.1:1, 2:5:20, 30:10:100, 200:100:1000, 2000:1000:10000]; 
maxUs=length(allUs);x=allUs; 
% f evaluated at each u point
fDoseResponse=zeros(1,maxUs);
for indU = 1:maxUs
    u=allUs(indU);
    fDoseResponse(indU)=doseResponse_MWC( LVal,n,vectCs,vectAlphas_i,u);
end 
% =======================================
% 2. Normalize curve
f_norm = fDoseResponse./fInf;
% =======================================
% 3. Fit to Hill Equation: x^n/(ec50^n+x^n)
F = @(c,x)(x.^c(1))./(c(2)^c(1)+x.^c(1)); % Hill function to fit to
% initial initial points for regression
x0=[0,1];  
% Numerical solver 
[c] = lsqcurvefit(F,x0,x,f_norm);
% c(1) = Hill number and c(2)=ec50; 
hillNo_pre=c(1);ec50 = c(2) ; 
% If the Hill number returned is imaginary or undefined, remove point by
% labeling as H=0
if isreal(hillNo_pre) == 0
    hillNo = 0 ;
else 
    hillNo = hillNo_pre;
end 

% Last update: 12/15/19 LL
end

