% This script generates and solves random instances of the 2-chain
% problem
%
% 6 Nov 2012
% J.Brooks

close all;
clear all;

OPTION = 5;
% 1 - new random data w/o delay
% 2 - new random data w/ delay
% 3 - new random N, same mus w/o delay (load file below)
% 4 - same N, mus, w/ random delay (load file below)
% 5 - same everything, but with M/M/1/K approximation (all others M/M/1)

if OPTION  > 2                          % load data
    load data/rand_2_6_20nov2012;
    
    if OPTION == 3
        entities = zeros(1,NUM_TESTS)*NaN;
    end

else                                    % else new data/structure
    type = [1 1 2 2 2 2 2 2];
    NUM_TESTS = 100;
    serviceRates = zeros(1,NUM_TESTS*length(type));
    entities = zeros(1,NUM_TESTS)*NaN;
end

cycles = find(type == 2);
chains = find(type == 1);

numCycles = length(cycles);
numChains = length(chains);
numStations = length(type);
decisionVars = numCycles*numChains;

% Output measures
cycleEquity = zeros(1,NUM_TESTS)*NaN;
centralEquity = zeros(1,NUM_TESTS)*NaN;
flowEquity = zeros(1,NUM_TESTS)*NaN;
centralDiff = zeros(1,NUM_TESTS)*NaN;
cycleDiff = zeros(1,NUM_TESTS)*NaN;
partition = zeros(1,NUM_TESTS)*NaN;
throughput = zeros(1,NUM_TESTS)*NaN;
centralCapacity = zeros(1,NUM_TESTS)*NaN;
cycleCapacity = zeros(1,NUM_TESTS)*NaN;

% Raw data
flows = zeros(1,NUM_TESTS*decisionVars);
travelDelay = zeros(1,NUM_TESTS*decisionVars);
waitTimes = zeros(1,NUM_TESTS*numStations);
utilizations = zeros(1,NUM_TESTS*numStations);
queueLengths = zeros(1,NUM_TESTS*numStations);

figure(1); subplot(311), title('wait times, 1p-2p'), hold on;
subplot(312), title('wait times, 2p-3p'), hold on;
subplot(313), title('wait times, 1p-3p'), hold on;

figure(2); subplot(311), title('wait times 1c-2c'), hold on;
subplot(312), title('wait times, 2c-3p'), hold on;
subplot(313), title('wait times, 1c-5p'), hold on;

for i = 1:NUM_TESTS
    if OPTION < 4
        N = randi(length(type)*2);      % number of entities
        entities(i) = N;
    else
        N = entities(i);                % get from file
    end
    
    if OPTION < 3
        mu = [randi(40, 1, length(chains))+5, randi(20,1,length(cycles)) ...
              + 10];                    % random service rates
    else                                % get from file
        mu = serviceRates((i-1)*numStations + 1 : i*numStations);
    end

    d = zeros(1,decisionVars); 
    
    if OPTION == 2 || OPTION == 4
        d = rand(1, decisionVars)*0.9+0.1; % U(0.1,1)
    end

    % Summary of random rates
    centralDiff(i) = max(mu(chains)) - min(mu(chains));
    cycleDiff(i) = max(mu(cycles)) - min(mu(cycles));

    % Solve
    if OPTION < 5
        [u,w,q,x,flag] = optimalAssignment(mu,type,d,N,1);
    else
        [u,w,q,x,flag] = optimalAssignment(mu,type,d,N,3);
    end

    
    % Test with proportional routing
    %    [u2,w2,q2,x2] = queueSim(N,mu,);

    % Calculate summary metrics
    cycleEquity(i) = max(u(cycles)) - min(u(cycles));
    flowEquity(i) = max(x) - min(x);
    centralEquity(i) = max(u(chains)) - min(u(chains));
    partition(i) = length(find(x < 0.01));
    throughput(i) = sum(x);
    centralCapacity(i) = throughput(i)/sum(mu(chains));
    cycleCapacity(i) = throughput(i)/sum(mu(cycles));
    
    % stash raw data
    flows((i-1)*decisionVars + 1 : i*decisionVars) = x;
    travelDelay((i-1)*decisionVars + 1 : i*decisionVars) = d;
    waitTimes((i-1)*numStations + 1 : i*numStations) = w;
    utilizations((i-1)*numStations + 1 : i*numStations) = u;
    queueLengths((i-1)*numStations + 1 : i*numStations) = q(1:end-1); ...
    % ignore travel
    
    serviceRates((i-1)*numStations + 1 : i*numStations) = mu;
    

    if ( flag == 1 ) & ( sum(x>0.01) == length(x) )%centralCapacity(i) < 0.9
        figure(1);        
        subplot(311), plot(sqrt(mu(cycles(1)))*w(cycles(1)), ...
                           sqrt(mu(cycles(2)))*w(cycles(2)), 'b.');
        subplot(312), plot(sqrt(mu(cycles(2)))*w(cycles(2)), ...
                           sqrt(mu(cycles(3)))*w(cycles(3)), 'b.');
        subplot(313), plot(sqrt(mu(cycles(1)))*w(cycles(1)), ...
                           sqrt(mu(cycles(3)))*w(cycles(3)), 'b.');
        figure(2);
        subplot(311), plot(sqrt(mu(chains(1)))*w(chains(1)), ...
                           sqrt(mu(chains(2)))*w(chains(2)), 'b.');
        subplot(312), plot(sqrt(mu(chains(2)))*w(chains(2)), ...
                           sqrt(mu(cycles(3)))*w(cycles(3)), 'b.');
        subplot(313), plot(sqrt(mu(chains(1)))*w(chains(1)), ...
                           sqrt(mu(cycles(5)))*w(cycles(5)), 'b.');
                       
%         figure(3);        
%         subplot(311), plot(sqrt(mu(cycles(1)))*w(cycles(1)), ...
%                            sqrt(mu(cycles(2)))*w(cycles(2)), 'b.');
%         subplot(312), plot(sqrt(mu(cycles(2)))*w(cycles(2)), ...
%                            sqrt(mu(cycles(3)))*w(cycles(3)), 'b.');
%         subplot(313), plot(sqrt(mu(cycles(1)))*w(cycles(1)), ...
%                            sqrt(mu(cycles(3)))*w(cycles(3)), 'b.');
    else
        figure(1);
        subplot(311), plot(sqrt(mu(cycles(1)))*w(cycles(1)), ...
                           sqrt(mu(cycles(2)))*w(cycles(2)), 'r.');
        subplot(312), plot(sqrt(mu(cycles(2)))*w(cycles(2)), ...
                           sqrt(mu(cycles(3)))*w(cycles(3)), 'r.');
        subplot(313), plot(sqrt(mu(cycles(1)))*w(cycles(1)), ...
                           sqrt(mu(cycles(3)))*w(cycles(3)), 'r.');
        figure(2);
        subplot(311), plot(sqrt(mu(chains(1)))*w(chains(1)), ...
                           sqrt(mu(chains(2)))*w(chains(2)), 'r.');
        subplot(312), plot(sqrt(mu(chains(2)))*w(chains(2)), ...
                           sqrt(mu(cycles(3)))*w(cycles(3)), 'r.');
        subplot(313), plot(sqrt(mu(chains(1)))*w(chains(1)), ...
                           sqrt(mu(cycles(5)))*w(cycles(5)), 'r.');
    end
end

figure, plot(cycleDiff, 'b.');
title(['Max. Cycle Service Rate ' ...
       'Difference']);

figure, plot(centralDiff, 'b.');
title(['Max. Central Service Rate ' ...
       'Difference']);

figure, plot(cycleEquity, 'b.');
title(['Max. Cycle Utilization Rate ' ...
       'Difference']);

figure, plot(centralEquity, 'b.');
title(['Max. Central Utilization Rate ' ...
       'Difference']);

figure, plot(flowEquity, 'b.');
title(['Max. Flow Rate ' ...
       'Difference']);

figure, plot(partition./length(x), 'b.');
title(['Percentage of flows ' ...
       'nominally zero']);

figure, plot(centralCapacity, 'b.'), hold on, 
plot(cycleCapacity, 'r.');
title('Capacities');

%const =   zeros(1,decisionVars); 
%for i = 1:length(chains)
%    for j = 1:length(cycles)
%        const((i-1)*length(cycles) + j) = mu(chains(i))*w(chains(i))^2 + mu(cycles(j))*w(cycles(j))^2 + d((i-1)*length(chains) + j)
%    end
%end
