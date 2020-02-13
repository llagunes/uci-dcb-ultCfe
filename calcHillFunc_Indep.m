function [ hillNo,xEC10,xEC90 ] = calcHillFunc_Indep(vs,d,k,n)
% This function calculates the Hill coefficient given the parameters above.
% This is for the Independent model

% Calculate fInf
fInf=(prod(vs))/(d+prod(vs));
% Calculate H
fun1 = @(u) doseResponse_Indep(d,vs,k,n,u)-0.1*fInf;% for EC10
fun2 = @(u) doseResponse_Indep(d,vs,k,n,u)-0.9*fInf;% for EC90
% fun3 = @(u) doseResponse_Indep(d,vs,k,n,u)-0.5*fInf;% for EC90

if n==1
    u1=[-0.1,10^10];u2=[0,10^10];
% elseif n==2
%     u1=[-1,10^3];u2=[0,10^10];
else 
     % u1=[-1,10^7];u2=[0,10^10];
     u1=[0,10^7];u2=[0,10^10];
end 


options=optimset('TolFun',1e-10);% tolerence
xEC10 = fzero(fun1,u1,options);
%xEC50= fzero(fun3,u2,options);
xEC90= fzero(fun2,u2,options);
hillNo = log(81)/log(xEC90/xEC10);

% For debugging
% doseResponse_MWC( LVal,n,vecEps,-10)-0.1*fInf
% doseResponse_MWC( LVal,n,vecEps,10)-0.1*fInf
% vs=vecEps;d=dVal;k=ks;

end


% % % if n>2
% % %     u1=[0, 100];
% % %     % u1=0;
% % %     u2 = [0,10000]; 
% % % elseif n==1
% % %     u1=[-1,10];u2=[0,10000];
% % % else
% % %     u1=[0,100];u2=[0,10^7];
% % % end



