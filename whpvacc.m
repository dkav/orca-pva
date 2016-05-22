% Orca population viability analysis using census counts
% Filename: whpvacc.m
% Author: Darren Kavanagh
% Date: 11.12.2008

clear all
load popds % file containing a list of datasets

% loop to fetch datasets and perform analysis
nods = (length(popds)); % no of datasets
results = zeros(nods,10);% matrix containing results for each dataset

for i=1:nods
 stwhpop = load(char(popds(i,1))); % loads as a structure
 fdnm = fieldnames(stwhpop); % get field names
 whpop = stwhpop.(fdnm{1});  % convert data from a structure to a matrix 
 clear stwhpop
 dataset = popds(i,2);  % dataset name

 curpop = (whpop(end,2)); % last population count 
 extinct = 10;            % extinction threshold

 interval = (length(whpop))-1; % the no of inter-census intervals
 duration = (whpop(end,1)) - (whpop(1,1)); % the duration of counts in yrs

 xdur = zeros(interval,1); % transformed duration variable
 ycnt = zeros(interval,1); % transformed count variable


 % Estimate µ and σ2 from a series of population censuses
 % ======================================================
 % Method outlined in Dennis et al. (1991)
 % parameter µ: determines how quickly the mean increases
 % parameter σ2: determines how quickly the variance in the normal
 %               distribution increases

 %  transform duration and count variables
 for j = 1:interval
     xdur(j) = sqrt(whpop(j+1,1)-whpop(j,1));
     ycnt(j) = (log(whpop(j+1,2)/whpop(j,2)))/xdur(j);
 end

 param1 = xdur\ycnt;    % slope of linear regression = parameter µ
 results(i,1) = param1; % store result
 pparam1 = [param1,0];  % set to polynomial as polyfit would return

 % plot transformed variables and regression line
 x = (0:0.2:max(xdur)+0.2); % calc x values for regression line
 y = (x * param1) + 0;      % calc y values for regression line
 figure (i)
     plot(xdur, ycnt, '*', x, y, '-')
     title (dataset)
     xlabel('x=√(t(j) - t(i))')
     ylabel('y=ln(N(j)/N(i))/√(t(j)-t(i))')
     box off

 ycntpred = polyval(pparam1,xdur); % calc predicated values
 resid = ycnt - ycntpred;          % calc residuals
 ssr = sum(resid.^2);              % calc squared sum of residuals
 param2 = ssr/(interval-1);        % mean square residual = parameter σ2
 results(i,2) = param2; % store result

 % Average Population Growth Rate
 % =============================
 cr = param1 + param2/2; % continuous rate of increase, r

 % lower 95% confidence limit for r
 lcr = cr + norminv(0.025) * sqrt(param2*((1/duration) ...
     + (param2/(2*(interval-1)))));
 % upper 95% confidence limit for r
 hcr = cr - norminv(0.025) * sqrt(param2*((1/duration) ...
     + (param2/(2*(interval-1)))));

 avgfr = exp(cr);   % average finite rate of increase, λ
 results(i,3) = avgfr; % store result
 avglfr = exp(lcr); % approximate lower 95% confidence limit for λ
 results(i,4) = avglfr; % store result
 avghfr = exp(hcr); % approximate upper 95% confidence limit for λ
 results(i,5) = avghfr; % store result


 % Probability of the Population Reaching a Lower Extinction Threshold
 % ===================================================================

 % slightly different estimate of parameter σ2
 nparam2 = (interval-1) * param2/interval;
 results(i,6) = nparam2; % store result

 % calc extinction probability
 if (param1<0)
     extprob =  1;
 else
     extprob = (extinct/curpop)^(2*param1/nparam2);
 end
 results(i,7) = extprob; % store result

 % The mean time to extinction
 % ===========================
 % Should only be used as a measure of extinction
 % risk when the estimate of µ is negative

 % log of the ratio between current populations and
 % the extinction threshold
 lncurext = log(curpop/extinct);

 if (param1<0)
     mtext = lncurext / (abs(param1));% mean time to extinction

     % lower 95% confidence limits
     if (mtext+norminv(0.025)*sqrt((lncurext^2)* ...
             nparam2/((param1^4)*duration))) < 0
         mtextl = 0;
     else
         mtextl = mtext+norminv(0.025)*sqrt((lncurext^2)* ...
             nparam2/((param1^4)*duration));
     end
     % upper 95% confidence limit
     mtexth = mtext-norminv(0.025)*sqrt((lncurext^2)* ...
         nparam2/((param1^4)*duration));
     
     results(i,8) = mtext;  % store result
     results(i,9) = mtextl;  % store result
     results(i,10) = mtexth; % store result
 end
 


 % The Cumulative Distribution Function (CDF) for
 % the Conditional Time to Extinction
 % ==================================================
 ftime = 5:5:1000; % series of times
 ftime = ftime';
 cdf = zeros(length(ftime),1);
 relviab = zeros(length(ftime),1);

 for j=1:(length(ftime))
     % function to calculate CDF
     cdf(j) = normcdf((-lncurext+abs(param1)*ftime(j)) ...
         /sqrt(nparam2*ftime(j))) ...
         + exp(2*lncurext*abs(param1)/nparam2)...
         * normcdf((-lncurext-abs(param1)*ftime(j)) ...
         /sqrt(nparam2*ftime(j)));
     relviab(j) = extprob * cdf(j); % calculate relative viability
 end

 relviab2 = zeros(length(ftime),1);
 cdf2 = zeros(length(ftime),1);
     cdf2 = normcdf((-lncurext+abs(param1)*ftime) ...
         /sqrt(nparam2*ftime)) ...
         + exp(2*lncurext*abs(param1)/nparam2)...
         * normcdf((-lncurext-abs(param1)*ftime) ...
         /sqrt(nparam2*ftime));
     relviab2 = extprob * cdf2; % calculate relative viability
 
 
% plot CDFs
 figure (nods+1)
    switch i
        case 1
            plot(ftime, relviab, 'b-')  
        case 2
            plot(ftime, relviab, 'r-')
        case 3
            plot(ftime, relviab, 'g-')   
    end
    hold on 
end

figure (nods+1)
    xlabel('Time (years)')
    ylabel('Cumulative Probability of Extinction')
    box off
    legend(char(popds(1,2)),char(popds(2,2)),char(popds(3,2)),'Location', 'SouthOutside', ...
             'Orientation','horizontal');
    legend('boxoff')
    axis([0 1000 0 1])
