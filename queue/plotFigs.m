close all;
clear all;

%--------------------
% Original, M/M/1
%--------------------
load data/rand_2_6_20nov2012.mat

cycles = find(type == 2);
chains = find(type == 1);

numCycles = length(cycles);
numChains = length(chains);
numStations = length(type);
numVars = numCycles*numChains;

maxThru = zeros(1,NUM_TESTS);
simThru = zeros(1,NUM_TESTS);
optThru = zeros(1,NUM_TESTS);

for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    idx2 = [(i-1)*numVars+1: i*numVars];
    maxThru(i) = min(sum(serviceRates(idx(cycles))), ...
                     sum(serviceRates(idx(chains))));
    optThru(i) = sum(flows(idx2));
end

figure(1); plot(entities, optThru./maxThru, 'b*');
hold on;

% Pure assignment
load data/rand_2_6_20nov2012_sim.mat

for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    simThru(i) = sum(throughputSIM(idx(chains)));
end
figure(1), plot(entities, simThru./maxThru, 'r^');
figure(2), plot(entities, simThru./optThru, 'r^');
hold on;

% Probabilistic assignment
load data/rand_2_6_20nov2012_mm1_sim_prob.mat

for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    simThru(i) = sum(throughputSIM(idx(chains)));
end
figure(1), plot(entities, simThru./maxThru, 'gs');
figure(2), plot(entities, simThru./optThru, 'gs');

%--------------------
% FEWER N
%--------------------
load data/rand_2_6_20nov2012_fewerN.mat
for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    idx2 = [(i-1)*numVars+1: i*numVars];
    maxThru(i) = min(sum(serviceRates(idx(cycles))), ...
                     sum(serviceRates(idx(chains))));
    optThru(i) = sum(flows(idx2));
end

figure(1); plot(entities, optThru./maxThru, 'b*');
hold on;

% Pure assignment
load data/rand_2_6_20nov2012_fewerN_mm1_sim.mat

for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    simThru(i) = sum(throughputSIM(idx(chains)));
end
figure(1), plot(entities, simThru./maxThru, 'r^');
figure(2), plot(entities, simThru./optThru, 'r^');
hold on;

% Probabilistic assignment
load data/rand_2_6_20nov2012_fewerN_mm1_sim_prob.mat

for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    simThru(i) = sum(throughputSIM(idx(chains)));
end
figure(1), plot(entities, simThru./maxThru, 'gs');
figure(2), plot(entities, simThru./optThru, 'gs');

figure(1);
legend({'Predicted', 'Pure', 'Prob.'}, 'location', 'southeast');
title('M/M/1 Approximation Throughput Results');
xlabel('Number of Entities');
ylabel('Relative Throughput');

figure(2);
legend({'Pure', 'Prob.'}, 'location', 'northeast');
title('M/M/1 Approximation Accuracy Results');
xlabel('Number of Entities');
ylabel('Actual Throughput Relative to Predicted');

clear all;

%--------------------
% Original, M/M/1/N
%--------------------
load data/rand_2_6_20nov2012_mm1k.mat

cycles = find(type == 2);
chains = find(type == 1);

numCycles = length(cycles);
numChains = length(chains);
numStations = length(type);
numVars = numCycles*numChains;

maxThru = zeros(1,NUM_TESTS);
simThru = zeros(1,NUM_TESTS);
optThru = zeros(1,NUM_TESTS);

for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    idx2 = [(i-1)*numVars+1: i*numVars];

    maxThru(i) = min(sum(serviceRates(idx(cycles))), ...
                     sum(serviceRates(idx(chains))));
    optThru(i) = sum(flows(idx2));
end

figure(3); plot(entities, optThru./maxThru, 'b*');
hold on;

% Pure assignment
load data/rand_2_6_20nov2012_mm1k_sim.mat

for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    simThru(i) = sum(throughputSIM(idx(chains)));
end
figure(3), plot(entities, simThru./maxThru, 'r^');
figure(4), plot(entities, simThru./optThru, 'r^');
hold on;

% Probabilistic assignment
load data/rand_2_6_20nov2012_mm1k_sim_prob.mat

for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    simThru(i) = sum(throughputSIM(idx(chains)));
end
figure(3), plot(entities, simThru./maxThru, 'gs');
figure(4), plot(entities, simThru./optThru, 'gs');

%--------------------
% FEWER N
%--------------------
load data/rand_2_6_20nov2012_fewerN_mm1k.mat
for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    idx2 = [(i-1)*numVars+1: i*numVars];
    maxThru(i) = min(sum(serviceRates(idx(cycles))), ...
                     sum(serviceRates(idx(chains))));
    optThru(i) = sum(flows(idx2));
end

figure(3); plot(entities, optThru./maxThru, 'b*');
hold on;

% Pure assignment
load data/rand_2_6_20nov2012_fewerN_mm1k_sim.mat

for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    simThru(i) = sum(throughputSIM(idx(chains)));
end
figure(3), plot(entities, simThru./maxThru, 'r^');
figure(4), plot(entities, simThru./optThru, 'r^');
hold on;

% Probabilistic assignment
load data/rand_2_6_20nov2012_fewerN_mm1k_sim_prob.mat

for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    simThru(i) = sum(throughputSIM(idx(chains)));
end
figure(3), plot(entities, simThru./maxThru, 'gs');
figure(4), plot(entities, simThru./optThru, 'gs');

figure(3);
legend({'Predicted', 'Pure', 'Prob.'}, 'location', 'southeast');
title('M/M/1/N Approximation Throughput Results');
xlabel('Number of Entities');
ylabel('Relative Throughput');

figure(4);
legend({'Pure', 'Prob.'}, 'location', 'northeast');
title('M/M/1/N Approximation Accuracy Results');
xlabel('Number of Entities');
ylabel('Actual Throughput Relative to Predicted');

clear all;

%--------------------
% Original, Finite-Source
%--------------------
load data/rand_2_6_20nov2012_finiteSource.mat

cycles = find(type == 2);
chains = find(type == 1);

numCycles = length(cycles);
numChains = length(chains);
numStations = length(type);
numVars = numCycles*numChains;

maxThru = zeros(1,NUM_TESTS);
simThru = zeros(1,NUM_TESTS);
optThru = zeros(1,NUM_TESTS);


for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    idx2 = [(i-1)*numVars+1: i*numVars];
    f = flows((i-1)*decisionVars + 1 : i*decisionVars);
    mus = serviceRates((i-1)*numStations + 1 : i*numStations);
    N = entities(i);
    maxThru(i) = min(sum(serviceRates(idx(cycles))), ...
                     sum(serviceRates(idx(chains))));
    for j=1:numChains
        centralFlow(j) = sum(f((j-1)*numCycles + 1:j* ...
                               numCycles));
        rho = centralFlow(j)/mus(chains(j));
        a = zeros(1,N);
        for k = 1:N
            a(k) = nchoosek(N,k)*factorial(k)*rho^k;
        end
        p0 = 1/(1+sum(a));
        delta = p0*sum([1:N].*a);

        centralL(j) = delta; %min(max(0,delta), N);
    end
    
    optThru(i) = sum(centralFlow.*(N-centralL));
end

figure(5); plot(entities, optThru./maxThru, 'b*');
hold on;

% Pure assignment
load data/rand_2_6_20nov2012_finiteSource_sim.mat

for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    simThru(i) = sum(throughputSIM(idx(chains)));
end
figure(5), plot(entities, simThru./maxThru, 'r^');
figure(6), plot(entities, simThru./optThru, 'r^');
hold on;

% Probabilistic assignment
load data/rand_2_6_20nov2012_finiteSource_sim_prob.mat

for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    simThru(i) = sum(throughputSIM(idx(chains)));
end
figure(5), plot(entities, simThru./maxThru, 'gs');
figure(6), plot(entities, simThru./optThru, 'gs');

%--------------------
% FEWER N
%--------------------
load data/rand_2_6_20nov2012_fewerN_finiteSource.mat
for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    idx2 = [(i-1)*numVars+1: i*numVars];
    f = flows((i-1)*decisionVars + 1 : i*decisionVars);
    mus = serviceRates((i-1)*numStations + 1 : i*numStations);
    N = entities(i);
    maxThru(i) = min(sum(serviceRates(idx(cycles))), ...
                     sum(serviceRates(idx(chains))));
    for j=1:numChains
        centralFlow(j) = sum(f((j-1)*numCycles + 1:j* ...
                               numCycles));
        rho = centralFlow(j)/mus(chains(j));
        a = zeros(1,N);
        for k = 1:N
            a(k) = nchoosek(N,k)*factorial(k)*rho^k;
        end
        p0 = 1/(1+sum(a));
        delta = p0*sum([1:N].*a);

        centralL(j) = delta; %min(max(0,delta), N);
    end
    
    optThru(i) = sum(centralFlow.*(N-centralL));
end

figure(5); plot(entities, optThru./maxThru, 'b*');
hold on;

% Pure assignment
load data/rand_2_6_20nov2012_fewerN_finiteSource_sim.mat

for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    simThru(i) = sum(throughputSIM(idx(chains)));
end
figure(5), plot(entities, simThru./maxThru, 'r^');
figure(6), plot(entities, simThru./optThru, 'r^');
hold on;

% Probabilistic assignment
load data/rand_2_6_20nov2012_fewerN_finiteSource_sim_prob.mat

for i = 1:NUM_TESTS
    idx = [(i-1)*numStations+1: i*numStations];
    simThru(i) = sum(throughputSIM(idx(chains)));
end
figure(5), plot(entities, simThru./maxThru, 'gs');
figure(6), plot(entities, simThru./optThru, 'gs');

figure(5);
legend({'Predicted', 'Pure', 'Prob.'}, 'location', 'southeast');
title('Finite-Source Approximation Throughput Results');
xlabel('Number of Entities');
ylabel('Relative Throughput');

figure(6);
legend({'Pure', 'Prob.'}, 'location', 'northeast');
title('Finite-Source Approximation Accuracy Results');
xlabel('Number of Entities');
ylabel('Actual Throughput Relative to Predicted');
