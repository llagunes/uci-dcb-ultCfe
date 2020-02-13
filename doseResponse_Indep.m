function [ f_indep,f_inf ] = doseResponse_Indep(d,vs,k,n,u)
% This function calculates the independent dose response given the
% parameters above

% Make I
pre_I = (dec2bin(0:(2^n)-1)=='1');[c1,b2]=sort(sum(pre_I,2));
I = pre_I(b2,:); % All cases of phosphoforms
n2 = size(I,1); % Number of phosphoforms to consider (2^n);

for jD=1:n2 
    % Define the activity level (Q_I) of protein 
    v_I = prod(vs(I(jD,:)==1));
    Q_I_diff(jD)=v_I/(d+v_I); % Activation function Q_diff = v^I/(d+v^I)
end

% Calculate Dose Response for different values of u
[ f_indep ] = calc_fU(Q_I_diff);
f_inf = prod(vs)/(d+prod(vs));% f_infinity
% =============================================
% Below is the functions calc_fU, which calculates f(u) at
% different u's: kinase concentrations
function [f_indep] = calc_fU(Q_I_diff)
    % This function calculates f(u), which is the Dose Response Curve
        % Q_I: vector containg protein activity.
        u_val = u;
        % First, determine S_I for n-sites
        p = u_val./(u_val+k); q = k./(k+u_val);
        for indR = 1:n2
            k_onTest = (I(indR,:)==1); k_offTest = I(indR,:)==0;
            prodOn = prod(p(k_onTest));  prodOff = prod(q(k_offTest)); %if isnan(prodOn); prodOn=0;end
            S_I(indR) = prodOn*prodOff; % Definition of S_I
            f_vals(indR) = S_I(indR)*Q_I_diff(indR); % Definition of f(u,v)
        end
        f_indep = sum(f_vals);


end

end 

