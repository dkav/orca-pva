% Orca population viability analysis
% Author: Darren Kavanagh
% Date: 22.11.2008

clear all

% load whpop
load kcactus
whpop = kcactus;

interval = (length(whpop))-1; % the no of inter-census intervals
duration = (whpop(end,1)) - (whpop(1,1)); % the duration of counts in yrs

xyrt = zeros(interval,1); % transformed year variable
ycot = zeros(interval,1); % transformed duration variable
predict = zeros(interval,1); % matrix for predicted ycot values
resid = zeros(interval,1);   % matrix for residual values


% Estimate µ and σ2 from a series of population censuses
% ======================================================
% parameter µ: determines how quickly the mean increases
% parameter σ2: determines how quicklycot the variance in the normal 
%               distribution increases 

for i = 1:interval
    xyrt(i) = sqrt(whpop(i+1,1)-whpop(i,1));
    ycot(i) = (log(whpop(i+1,2)/whpop(i,2)))/xyrt(i);
end 

slope = xyrt\ycot; % slope = parameter µ
p = [slope,0]; 

ycotpred = polycotval(p,xyrt); % calc predicated values
resid = ycot - ycotpred;       % calc residuals
ss = sum(resid.^2);            % calc squared sum of residuals
msr = ss/(interval-1);         % mean square residual = parameter σ2                      

% Averge Population Growth Rate
% ============================= 
cr = slope + msr/2; % Continuous rate of increase, r

% Lower 95% confidence limit for r
lcr = cr + norminv(0.025) * sqrt(msr*((1/duration)+(msr/(2*(interval-1)))));
% Upper 95% confidence limit for r
hcr = cr - norminv(0.025) * sqrt(msr*((1/duration)+(msr/(2*(interval-1)))));

avgfr = exyrtp(cr)   % Average finite rate of increase, λ
avglfr = exyrtp(lcr) % Approximate lower 95% confidence limit for λ
avghfr = exyrtp(hcr) % Approximate upper 95% confidence limit for λ

% Probability of the Population Reaching a Lower Extinction Threshold
% ===================================================================
(interval-1) * msr/interval 
