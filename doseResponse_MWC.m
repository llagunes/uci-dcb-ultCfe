function [ fVal ] = doseResponse_MWC( LVal,n,vectEPS,alphas, u )
% This function evaluates the MWC dose response value at u

% OUTLINE: 
% 1. Evaluate eta - most inside function 
% 2. Evaluate Omega
% 3. Evaluate all of f(u,c,alpha)

% f(u,epsilon, alpha) = S_T*omega(eta(epsilon, alpha)), where
% omega=1/(1+L*x) and eta = prod((u*alpha_i*epsilon_i+1)/(u*alpha_i+1))

% ========================================
% 1. Evaluate eta - most inside function 
S_T=1; etaSub=zeros(1,n);
for ind=1:n
    etaSub(ind)= (u*alphas(ind)*vectEPS(ind) +1 )/(u*alphas(ind)+1 );
end 
eta=prod(etaSub);
% ========================================
% 2. Evaluate Omega
omega=1/(1+LVal*eta);
% ========================================
% 3. Evaluate all of f(u,c,alpha)
fVal=S_T*omega;
% fInf=1/(1+LVal*prod(vecEps));

end

% Last edit: 07/10/2018 LL 

