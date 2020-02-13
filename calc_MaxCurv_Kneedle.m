function [ x_locKnee, y_knee] = calc_MaxCurv_Kneedle( xVals, yVals,cutOffval,mtdi,top)
% This function will locate the "knee" of a curve.

% mtdi=1
% by finding the location where the curve has the same slope as the line that
% corosses the first and last point (defined by middle percentages) 

% mtdi=2
% by finding the location where the curvature is minimum

% mtdi=3
% by following paper: Finding a Kneedle in a Haystack: Detecting Knee Points in System Behavior
% link: https://raghavan.usc.edu/papers/kneedle-simplex11.pdf
% Specifically only Figure 2a 

% ======================================
% Outline: 
% 1. Read in curve based on cut of values (the middle cutOff% of points) 
% 2. Calculate the location of the knee based on which method works best
% (mtdi) 
% 3. Export the knees

% =======================================
% 1. Read in curve based on cut off values (the middle cutOff% of points)

% Take the middle cutOffval% of data points
totPts=size(xVals,1); 

% if top percent of data is to be used, need to take top instead of middle
if top == 0 % use middle 
    ind1stPt = ceil((1-cutOffval)/2*totPts);indLstPt = totPts-ind1stPt;
    % define points
    x=xVals(ind1stPt+1:indLstPt);
    y=yVals(ind1stPt+1:indLstPt);
elseif top == 1 % use top percent 
    ind1stPt = ceil((1-cutOffval)/2*totPts);
    % flip order since in descending order right now
    [xpre, indxs]=sort(xVals,'ascend');[ypre, indys]=sort(yVals,'ascend');
    x=xpre(ind1stPt+1:end);
    y=ypre(ind1stPt+1:end);
elseif top == 2 % use bottom percent
    ind1stPt = ceil((1-cutOffval)/2*totPts);indLstPt = totPts-ind1stPt;
    x=xVals(1:indLstPt);
    y=yVals(1:indLstPt);
end 


% ======================================
% 2. Calculate the location of the knee based on which method works best
if mtdi==1 % slopes
    % Find line between first and last point of graph
    x1=x(1); x2=x(end); y1=y(1); y2=y(end);
    m =(y2-y1)/(x2-x1); % slope of line
    dyNum=diff(y)./diff(x); % Find slope of curve
     % Find where y' = m
     diffSlopes=abs(dyNum-m); [knee,xKnee]=min(diffSlopes);
     % Knee where slopes are closest
     xp1=x(xKnee+1);yp1=y(xKnee+1);
elseif mtdi==2 % curvature
    % Find y'
    dyNum=diff(y)./diff(x);
    % Find y'' using central differences
    gradX=diff(x);A=y';N=size(A,2);
    for j=2:N-1
        G(:,j) = ( A(:,j+1) -2*A(:,j) + A(:,j-1) );
    end
    y2p=(G')./(gradX.^2);
    kappa = abs(y2p./(1+dyNum.^2).^(3/2));% calculate curvature kappa
    % Knee where kappa is smallest
    [knee,xKnee]=max(kappa); xp1=x(xKnee);yp1=y(xKnee);
elseif mtdi==3 % Sapotta
    % Find line between first and last point of graph
    x1=x(1); x2=x(end); y1=y(1); y2=y(end);
    m= (y2-y1)/(x2-x1); b=y1-m*x1;line=m*x+b;
    % Calculate the difference between each data point and the line
    diffPts2Line=abs(y - line);
    % Knee is where the max difference exists
    [knee,xKnee]=max(diffPts2Line);
    xp1=x(xKnee);yp1=y(xKnee);   
end 
% ====================================
% 3. Export the knees
x_locKnee =  xp1;
y_knee = yp1;

% =======================================



% Last update: 01/20/19 LL
end

