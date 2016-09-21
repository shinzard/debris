% This is the calling script to run the simulation experiments. 
% 
% 4 Nov 2013
% J.Brooks
% 
% Last Modified: 4 Nov 2013
% 
clear all;
close all;

addpath('../');                         % to get optimization functions

% Parameters
dispatchParams.method = 'Bounded';      % {'Fixed', 'Bounded', 'Optimal'}
dispatchParams.inter = 0.44;            % hours
dispatchParams.forgetFactor = 0.95;     % EWMA forgetting factor
dispatchParams.Update = false;          % update prices optimally?

compParams.deltaA = -0.3;               % composition effect
compParams.deltaN = -2.2;               % composition effect
compParams.distRW = true;               % distances random walk?

% DoE
reps = 20;                              % number of replications
duration = 12*60;                       % length of simulation (hours)

% Network topology
mus = [10.26, 6.62, 6, 6, 6, 6, 6, 6];  % for parallel servers,
                                        % this is 'base mu'
type = [1 1 2 2 2 2 2 2];               % NOTE: simulateQ assumes
                                        % that the central servers
                                        % are listed first (type=1)
% Initial travel times (hrs)
dij = [ 0.17 0.62; ...                  % from prices.xlsx       
        0.23 0.61; ...
        0.24 3.35; ...
        0.22 0.76; ...
        0.47 0.49; ...
        0.50 1.11 ];

% Remote --> Central routing probabilities (Debris type) ....these
% could also be random if desired..
pr = [ 0.33 0.67; ...
       0.60 0.40; ...
       0.92 0.08; ...
       0.20 0.80; ...
       0.92 0.08; ...
       0.43 0.57 ];

% Pricing strategies
 % priceOpt = [1 1 1 1 1 1; ...              % uniform
 %             1, 0.9345, 0.9515, 1.0875, 0.9438, 1.1404;... % optimal
 %             0.85, 0.69, 0.91, 1.18, 0.85, 1.53; ...% time
 %             1.70, 0.74, 0.29, 1.93, 0.28, 1.06];% distance
 
% On the simplex
priceOpt = [ 0.1667    0.1667    0.1667    0.1667    0.1667    0.1667; % uniform
%            0.1651    0.1543    0.1571    0.1795    0.1558    0.1883; % optimal
             0.1414    0.1148    0.1514    0.1963    0.1414    0.2546; % time
             0.2833    0.1233    0.0483    0.3217    0.0467    0.1767]; % distance

numPrices = size(priceOpt,1);

% Starting allocation (uniform)
numTrucks = 24;
numDrts = sum(type==2);
N = repmat(numTrucks/numDrts, 1, numDrts);

run = 1;
throughput = zeros(1,2*numPrices*reps + 3*reps);
revenue = zeros(1,2*numPrices*reps + 3*reps);
inequity = zeros(1,2*numPrices*reps + 3*reps);
inequity2 = zeros(1,2*numPrices*reps + 3*reps);
inequity3 = zeros(1,2*numPrices*reps + 3*reps);
inequity4 = zeros(1,2*numPrices*reps + 3*reps);
T = zeros(1,2*numPrices*reps + 3*reps);
meanFam = zeros(1,2*numPrices*reps + 3*reps);
meanMu = zeros(1,2*numPrices*reps + 3*reps);
maxSize = zeros(1,2*numPrices*reps + 3*reps);
priceStrategy = zeros(1,2*numPrices*reps + 3*reps);

for temp = [0.73, 0.2]
    dispatchParams.T = temp            % level of rationality

    for p = 1:numPrices
        prices = priceOpt(p,:)         % prices
    
        for i = 1:reps
            fprintf('.');
            [u,w,q,x,m,s,f] =  simulateQ(N, mus, type, -1, dij, dispatchParams, ...
                                   compParams, pr, duration, 0, prices);
            
            throughput(run) = sum(x(find(type==1)));
            revenue(run) = sum(prices.*x(type==2));
            inequity(run) =  sqrt(sum(((x(type==2)./mean(x(type==2))-1).^2))/(length(x)*length(x)-1));% Flow
            inequity2(run) = sqrt(sum(((w(type==2)./mean(w(type==2))-1).^2))/(length(x)*length(x)-1));% Wait
            inequity3(run) = sqrt(sum(((q(type==2)./mean(q(type==2))-1).^2))/(length(x)*length(x)-1));% Length
            inequity4(run) = sqrt(sum(((u(type==2)./mean(u(type==2))-1).^2))/(length(x)*length(x)-1));% Utilization
            T(run) = temp;
            meanFam(run) = mean(f);
            meanMu(run) = mean(m);
            maxSize(run) = max(s);
            priceStrategy(run) = p;
            run = run + 1;
        end
    end
end

% % Just run optimal prices...
prices = [1, 0.9345, 0.9515, 1.0875, 0.9438, 1.1404];% initial optimal
prices = prices/sum(prices);            % (on simplex)
dispatchParams.Update = true;          % but, update prices

for temp = [0.73, 0.2]
    dispatchParams.T = temp            % level of rationality
    for i = 1:reps
        fprintf('.');
        [u,w,q,x,m,s,f] =  simulateQ(N, mus, type, -1, dij, dispatchParams, ...
                               compParams, pr, duration, 0, prices);
        
        throughput(run) = sum(x(find(type==1)));
        revenue(run) = sum(prices.*x(type==2));
        inequity(run) =  sqrt(sum(((x(type==2)./mean(x(type==2))-1).^2))/(length(x)*length(x)-1));% Flow
        inequity2(run) = sqrt(sum(((w(type==2)./mean(w(type==2))-1).^2))/(length(x)*length(x)-1));% Wait
        inequity3(run) = sqrt(sum(((q(type==2)./mean(q(type==2))-1).^2))/(length(x)*length(x)-1));% Length
        inequity4(run) = sqrt(sum(((u(type==2)./mean(u(type==2))-1).^2))/(length(x)*length(x)-1));% Utilization
        T(run) = temp;
        meanFam(run) = mean(f);
        meanMu(run) = mean(m);
        maxSize(run) = max(s);
        priceStrategy(run) = numPrices+1;
        run = run + 1;
    end
end

% Now just fixed optimal partition...
 dispatchParams.method = 'Optimal'      % {'Fixed', 'Bounded',
                                        % 'Optimal'}
 dispatchParams.Update = false;         % just to save computation

 for i = 1:reps
     fprintf('.');
     [u,w,q,x,m,s,f] =  simulateQ(N, mus, type, -1, dij, dispatchParams, ...
                            compParams, pr, duration, 0, prices);
    
     throughput(run) = sum(x(find(type==1)));
     revenue(run) = sum(prices.*x(type==2));
     inequity(run) =  sqrt(sum(((x(type==2)./mean(x(type==2))-1).^2))/(length(x)*length(x)-1));% Flow
     inequity2(run) = sqrt(sum(((w(type==2)./mean(w(type==2))-1).^2))/(length(x)*length(x)-1));% Wait
     inequity3(run) = sqrt(sum(((q(type==2)./mean(q(type==2))-1).^2))/(length(x)*length(x)-1));% Length
     inequity4(run) = sqrt(sum(((u(type==2)./mean(u(type==2))-1).^2))/(length(x)*length(x)-1));% Utilization
     T(run) = NaN;
     meanFam(run) = mean(f);
     meanMu(run) = mean(m);
     maxSize(run) = max(s);
     priceStrategy(run) = numPrices+2;
     run = run + 1;
 end

