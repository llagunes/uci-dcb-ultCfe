function [ actualEpsVec,meanEps,allCvs,hillsOrd,hillsOrdSame,vectAllKs ] = hillsVsMeanCV_INDEP( n, maxRuns, rangV,rangKs,dVal )
% This function will return the mean(v), cv(v), and Hill
% Number based on the ranges available for v

% Ranges of v and k
vMaxDeg=rangV(2); vMinDeg=rangV(1);
KMaxDeg=rangKs(2);KMinDeg=rangKs(1);
% Store all the v's, and Hill numbers
varVs=zeros(maxRuns,n);vMeans=zeros(maxRuns,1);allCvs_pre=zeros(maxRuns,1);
allHs=zeros(maxRuns);allHs2=zeros(maxRuns);
% make vector of k values for this n 
rValpha = KMinDeg + (KMaxDeg-KMinDeg).*rand(maxRuns,n); % randomly chose the degree's
vectAllKs=10.^rValpha; % alpha - modification rate parameter
% Determine how H is affected by mean(epsilon)
for v1=1:maxRuns % For each v vector
    r = vMinDeg + (vMaxDeg-vMinDeg).*rand(n,1);% randomly sample epsilon
    vectEPS=10.^r;
    vectKs=vectAllKs(v1,:);
    % Calculate mean of this vector and make 2nd vector of v_i=mean
    newMean=mean(vectEPS);
    vectEps2=newMean*ones(n,1);
    % Store epsilon vector, mean, and cv
    varVs(v1,:)=vectEPS;vMeans(v1)=newMean;
    allCvs_pre(v1)=std(vectEPS)/mean(vectEPS);
    % Calculate H for epsilon1 and epsilon2
    hillNo = calcHillFunc_Indep(vectEPS,dVal,vectKs,n);
    hillNoSame = calcHillFunc_Indep(vectEps2,dVal,vectKs,n);
    % Each column is H for each trial, row is H for each combo of epsilons
    allHs(v1)=hillNo;
    allHs2(v1)=hillNoSame;
end

% Sort in ascending mean order
[listed,indOrd]=sort(vMeans,'ascend');
meanEps=listed; hillsOrd=allHs(indOrd);hillsOrdSame=allHs2(indOrd);
allCvs=allCvs_pre(indOrd);

actualEpsVec=varVs;




% Last edit: 02/13/20 LL


end

