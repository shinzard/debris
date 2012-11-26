% Uncomment for execution on morse
%addpath('/usr/share/octave/packages/3.2/queueing-1.1.1');
close all;
clear all;
more off;

% These need to match the generated data from optimalAssignment
%load rand_2_6_20nov2012
%load rand_2_6_20nov2012_sim

load rand_2_6_20nov2012_fewerN
load rand_2_6_20nov2012_sim_fewerN

%type = [1 1 2 2 2 2 2 2]; IN MAT FILE
%NUM_TESTS = 100; IN MAT FILE

cycles = find(type == 2);
chains = find(type == 1);

numCycles = length(cycles);
numChains = length(chains);
numStations = length(type);
decisionVars = numCycles*numChains;

%throughputSIM = zeros(1,NUM_TESTS*decisionVars);
%waitTimesSIM = zeros(1,NUM_TESTS*numStations);
%utilizationsSIM = zeros(1,NUM_TESTS*numStations);
%lengthsSIM = zeros(1,NUM_TESTS*numStations);
    %     class, -1 for rational routing)

theoreticalMax = zeros(1,NUM_TESTS);
optimalMax = zeros(1,NUM_TESTS);
simThroughput = zeros(1,NUM_TESTS);

for i = 1:NUM_TESTS
    f = flows((i-1)*decisionVars + 1 : i*decisionVars);
    d = travelDelay((i-1)*decisionVars + 1 : i*decisionVars);
    mus = serviceRates((i-1)*numStations + 1 : i*numStations);
    xt = throughputSIM((i-1)*numStations + 1 : i*numStations);
    
    theoreticalMax(i) = min(sum(mus(chains)), sum(mus(cycles)));
    optimalMax(i) = sum(f);
    simThroughput(i) = sum(sum(xt(:,chains)));
end

figure, plot(entities, optimalMax./theoreticalMax, 'b*', 'markersize', ...
             10);
hold on;
plot(entities, simThroughput./theoreticalMax, 'r^', 'markersize', ...
     8);
title('Throughput relative to theoretical maximum');
legend({'Approx Optimization', 'Simulation w/ Pure Assignment'}, ...
       'location', 'southeast');

figure, plot(entities, simThroughput./optimalMax, 'b*', 'markersize', ...
             10);
title('Relative simulation throughput relative to approx optimal');