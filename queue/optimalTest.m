% This script generates and solves random instances of the 2-chain
% problem
%
% 6 Nov 2012
% J.Brooks

close all;
clear all;

type = [1 1 1 2 2 2 2 2 2 2 2 2 2 2];
cycles = find(type == 2);
chains = find(type == 1);
d = 0;
NUM_TESTS = 100;

% Output measures
cycleEquity = zeros(1,NUM_TESTS)*NaN;
centralEquity = zeros(1,NUM_TESTS)*NaN;
flowEquity = zeros(1,NUM_TESTS)*NaN;
centralDiff = zeros(1,NUM_TESTS)*NaN;
cycleDiff = zeros(1,NUM_TESTS)*NaN;
partition = zeros(1,NUM_TESTS)*NaN;
entities = zeros(1,NUM_TESTS)*NaN;
throughput = zeros(1,NUM_TESTS)*NaN;
centralCapacity = zeros(1,NUM_TESTS)*NaN;
cycleCapacity = zeros(1,NUM_TESTS)*NaN;

figure(1); subplot(311), title('wait times, 1-2'), hold on;
subplot(312), title('wait times, 2-3'), hold on;
subplot(313), title('wait times, 1-3'), hold on;

for i = 1:NUM_TESTS
    N = randi(length(type)*5);          % number of entities
    entities(i) = N;
    mu = randi(20, 1, length(type))+5;  % 5-25 service rates
    centralDiff(i) = max(mu(chains)) - min(mu(chains));
    cycleDiff(i) = max(mu(cycles)) - min(mu(cycles));
    [u,w,q,x,lagrange] = optimalAssignment(mu,type,d,N);
    cycleEquity(i) = max(u(cycles)) - min(u(cycles));
    flowEquity(i) = max(x) - min(x);
    centralEquity(i) = max(u(chains)) - min(u(chains));
    partition(i) = length(find(x < 0.1));
    throughput(i) = sum(x);
    centralCapacity(i) = throughput(i)/sum(mu(chains));
    cycleCapacity(i) = throughput(i)/sum(mu(cycles));
    figure(1);
    if centralCapacity(i) < 0.9
        subplot(311), plot(sqrt(mu(cycles(1)))*w(cycles(1)), sqrt(mu(cycles(2)))*w(cycles(2)), 'b.');
        subplot(312), plot(sqrt(mu(cycles(2)))*w(cycles(2)), sqrt(mu(cycles(3)))*w(cycles(3)), 'b.');
        subplot(313), plot(sqrt(mu(cycles(1)))*w(cycles(1)), sqrt(mu(cycles(3)))*w(cycles(3)), 'b.');
    else
        subplot(311), plot(sqrt(mu(cycles(1)))*w(cycles(1)), sqrt(mu(cycles(2)))*w(cycles(2)), 'r.');
        subplot(312), plot(sqrt(mu(cycles(2)))*w(cycles(2)), sqrt(mu(cycles(3)))*w(cycles(3)), 'r.');
        subplot(313), plot(sqrt(mu(cycles(1)))*w(cycles(1)), sqrt(mu(cycles(3)))*w(cycles(3)), 'r.');
    end
end

figure, plot(cycleDiff, 'b.'), title(['Max. Cycle Service Rate ' ...
                    'Difference']);
figure, plot(centralDiff, 'b.'), title(['Max. Central Service Rate ' ...
                    'Difference']);
figure, plot(cycleEquity, 'b.'), title(['Max. Cycle Utilization Rate ' ...
                    'Difference']);
figure, plot(centralEquity, 'b.'), title(['Max. Central Utilization Rate ' ...
                    'Difference']);
figure, plot(partition./length(x), 'b.'), title(['Percentage of flows ' ...
                    'nominally zero']);

figure, plot(centralCapacity, 'b.'), hold on, 
plot(cycleCapacity, 'r.');
title('Capacities');
