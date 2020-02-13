function [ hillNo,xEC50 ] = calcHillFunc_Levitzki(vectCs,LVal,n,vectAlphas_i)
% This function will return the hill number according to the Levitsky
% definition in 1978 

% Outline: 
% 1. Find the EC50 - enzyme concentration at which there is a 50% response
% 2. Calculate the slope around the EC50 
% 3. Calculate n50 = 4*EC50*f'(EC50)

% =============================================
% For numerically finding EC50
% Calculate fInf and f for fzero solver
fInf=1/(1+LVal*prod(vectCs));
fun3 = @(u)( doseResponse_MWC( LVal,n,vectCs,vectAlphas_i,u)-0.5*fInf);% for EC90
   
% =============================================
% 1. Find the EC50 - enzyme concentration at which there is a 50% response
u2=[0,10^7];% initial point for solver
options=optimset('TolFun',1e-10);% tolerence
% Test points
tp1=doseResponse_MWC( LVal,n,vectCs,vectAlphas_i,0)-0.5*fInf;
tp2=doseResponse_MWC( LVal,n,vectCs,vectAlphas_i,10^7)-0.5*fInf;

if sign(tp1)==sign(tp2)
    u2=[-1,10^7];
else 
    u2=u2;
end 
xEC50= fzero(fun3,u2,options);
% =============================================
% 2. Calculate the slope around the EC50 
dp=0.0001; 
p2=doseResponse_MWC( LVal,n,vectCs,vectAlphas_i,xEC50)/fInf;
p1=doseResponse_MWC( LVal,n,vectCs,vectAlphas_i,(xEC50+dp))/fInf;
m=(p1-p2)/dp;
% =============================================
% 3. Calculate n50 = 4*EC50*f'(EC50)
n_50=4*xEC50*m;
% If the Hill number returned is imaginary or undefined, remove point by
% labeling as H=0
if n_50 < 0
    hillNo = 0 ;
else 
    hillNo = n_50;
end 

% Last update: 12/15/19 LL
end

