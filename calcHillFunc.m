function [ hillNo,xEC10,xEC90 ] = calcHillFunc( vectCs,LVal,n,vectAlphas_i )
% This function calculates the Hill coefficient given the parameters above.
% Using the Goldbeter-Koshland definition for H = ln(81)/ln(EC90/EC10)

% Outline: 
% 1. Find EC10 and EC90 - enzyme concentration at which there is a 10, 90% response
% 2. Calculate H = ln(81)/ln(EC90/EC10)

% For numerically finding EC10,90
% Calculate fInf
fInf=1/(1+LVal*prod(vectCs));
% Calculate f for fzero solver
fun1 = @(u)( doseResponse_MWC( LVal,n,vectCs,vectAlphas_i,u)-0.1*fInf);% for EC10
fun2 = @(u)( doseResponse_MWC( LVal,n,vectCs,vectAlphas_i,u)-0.9*fInf);% for EC90
% Test points for sover
if n==1
    u1=[-1,10^10];u2=[-1,10^10];
elseif n==2
    testPtU=0; testPtU2=-2;
    yTest=doseResponse_MWC( LVal,n,vectCs,vectAlphas_i,testPtU)-0.1*fInf; 
    yTest2=doseResponse_MWC( LVal,n,vectCs,vectAlphas_i,testPtU2)-0.1*fInf;    
    if yTest<0
        u1=[0,10^7];u2=[0,10^10];
    elseif yTest2 < 0
        u1=[-2,10^7];u2=[0,10^10];
    elseif (vectCs(1)>=1 && vectCs(2) >=0.5) || (vectCs(2)>=1 && vectCs(1) >=0.5)
            u1=[-1.5,10^7];u2=[0,10^10]; 
    elseif (vectCs(1)>=0.5 && vectCs(2) >= 0.5) || (vectCs(2)>=0.5 && vectCs(1) >= 0.5)
        u1=0;u2=1; 
    else
        u1=0;u2=[0,10^10];
    end 
else
    testPtU=0; testPtU2=-1;
    yTest=doseResponse_MWC( LVal,n,vectCs,vectAlphas_i,testPtU)-0.1*fInf; 
    yTest2=doseResponse_MWC( LVal,n,vectCs,vectAlphas_i,testPtU2)-0.1*fInf;    
    if yTest<0
        u1=[0,10^7];u2=[0,10^10];
    elseif yTest2 < 0
        u1=[-1,10^7];u2=[-1,10^10];
    else
        u1=0;u2=[0,10^10];
    end 
end 

% 1. Find EC10 and EC90 - enzyme concentration at which there is a 10, 90% response
options=optimset('TolFun',1e-10);% tolerence
xEC10 = fzero(fun1,u1,options);
xEC90= fzero(fun2,u2,options);
% 2. Calculate H = ln(81)/ln(EC90/EC10)
% If the Hill number returned is imaginary or undefined, remove point by
% labeling as H=0
if (xEC10 <= 0 ) || (xEC90 <=0)
    hillNo = 0 ;
else 
    hillNo = log(81)/log(xEC90/xEC10);
end 

% Last update: 12/12/19 LL
end

