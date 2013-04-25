% This script generates and solves random instances of the 2-chain
% problem
%
% 6 Nov 2012
% J.Brooks

close all;
clear all;

OPTION = 2;
PLOT = false;
% 1 - new random data w/o delay
% 2 - new random data w/ delay
% 3 - new random N, same mus w/o delay (load file below)
% 4 - same N, mus, w/ random delay (load file below)
% 5 - same everything, but with new approximation method (all others M/M/1)

if OPTION  > 2                          % load data
    load data/rand_2_6_20nov2012_fewerN;
    
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
exitFlags = zeros(1,NUM_TESTS)*NaN;
% EQUITY
perfDecrease = zeros(1,NUM_TESTS)*NaN;
equityIncrease = zeros(1,NUM_TESTS)*NaN;
termFlag = zeros(1,NUM_TESTS)*NaN;

% Raw data
flows = zeros(1,NUM_TESTS*decisionVars);
travelDelay = zeros(1,NUM_TESTS*decisionVars);
waitTimes = zeros(1,NUM_TESTS*numStations);
utilizations = zeros(1,NUM_TESTS*numStations);
queueLengths = zeros(1,NUM_TESTS*numStations);

if PLOT
    figure(1); subplot(311), title('wait times, 1p-2p'), hold on;
    subplot(312), title('wait times, 2p-3p'), hold on;
    subplot(313), title('wait times, 1p-3p'), hold on;

    figure(2); subplot(311), title('wait times 1c-2c'), hold on;
    subplot(312), title('wait times, 2c-3p'), hold on;
    subplot(313), title('wait times, 1c-5p'), hold on;
end

for i = 1:NUM_TESTS
    disp(sprintf('i=%d',i));
    if OPTION < 4
        N = randi(length(type)*5);      % number of entities
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
        
        weights = repmat(1,1,decisionVars);
        lastW = zeros(1,decisionVars);
        lastX = zeros(1,decisionVars);
        for j = 1:100
            [u,w,q,x,f] = optimalAssignmentWeights(mu,type,d,weights,N, ...
                                               1);
            inequity = 1/decisionVars*sqrt((x./mean(x)-1).^2);
            if j == 1
                initPerf = sum(x);
                initI = sum(inequity);
            end
            neg = (x<mean(x)) - (x>mean(x));
            weights = weights+inequity.*neg*1/(decisionVars);
            if PLOT
                figure(50), subplot(211), plot([1:decisionVars]*2+0.01*(j-1), ...
                                               x, 'b.'), hold on;
                subplot(212), plot([1:decisionVars]*2+0.01*(j-1), weights, ...
                                   'r.'), hold on;
            
                figure(51), subplot(311), plot(sum(inequity),sum(x), ...
                                               'b.'), hold on;
                subplot(312), plot(j, sum(x), 'b.'), hold on;
                subplot(313), plot(j, sum(inequity),'r.'), hold on;
            end
            
            if norm(weights-lastW) < 0.02
                flag = 1;
                break;
            end
            
            if  norm(x-lastX) < 0.05
                flag = 2;
                break;
            end
            lastW = weights;
            lastX = x;
            flag = 0;
        end
        [Nvec, p] = routing(type, mu, d, x, N, 1);
        if PLOT
            figure(100), plot(Nvec,x./Nvec, 'b.'), hold on;
            title('Individual flow by team size');
            figure(101), plot(mu,x./Nvec, 'b.'), hold on;
            title('Individual flow by loading rate');
        end
    else
        [u,w,q,x,f] = optimalAssignment(mu,type,d,N,3);
    end

    perfDecrease(i) = sum(x)/initPerf;
    equityIncrease(i) = sum(inequity)/initI;
    termFlag(i) = flag;
    
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
    queueLengths((i-1)*numStations + 1 : i*numStations) = q(1:end-1);
    exitFlags(i) = f;
    
    % ignore travel
    
    serviceRates((i-1)*numStations + 1 : i*numStations) = mu;
    

    if ( flag == 1 ) & ( sum(x>0.01) == length(x) ) & PLOT%centralCapacity(i) < 0.9
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
    elseif PLOT
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

if PLOT
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
end

%const =   zeros(1,decisionVars); 
%for i = 1:length(chains)
%    for j = 1:length(cycles)
%        const((i-1)*length(cycles) + j) = mu(chains(i))*w(chains(i))^2 + mu(cycles(j))*w(cycles(j))^2 + d((i-1)*length(chains) + j)
%    end
%end
