% Stochastic Population Model
% Filename: stochpop.m
% Author: Darren Kavanagh
% Date: 11.12.2008

clear all

stpop = 10;             % starting population
tinter = 21;            % time intervals
maxtime = tinter-1;     % final time
timespan = 0:1:maxtime; % timespan vector
avgr = 1.05;            % average growth rate
varr = 0.015;           % variance of average growth rate
nosim = 20;             % no of simulations

% matrix to store all simulation results
allN = zeros(tinter, nosim); 

% average population loop
N = zeros(tinter,1);
N(1)= stpop;  
for t=1:maxtime        
        N(t+1) = avgr * N(t);             % calc new population
end
avgN = N

% simulation loop
for i=1:nosim
    N = zeros(tinter,1);
    N(1)= stpop;            % set starting population
    % time interval loop
    for t=1:maxtime
        r = avgr + sqrt(varr) * randn; % get a random growth rate
        N(t+1) = r * N(t);             % calc new population
    end
    allN(:,i) = N;
end

% plot results
figure (1)
    box off    
    plot(timespan, avgN, 'r-*',...
         timespan, max(allN,[],2),'b-s',...
         timespan, min(allN,[],2), 'b-o', ...
     timespan,allN', 'b-')
    xlabel('Time, t');
    ylabel('Populations Size, N(t)');
    legend('Mean λ', 'Max λ', 'Min λ','Location', 'SouthOutside', ...
           'Orientation','horizontal');
    legend('boxoff')