% This script is intended to determine the approximate delay
% function (and flow) for a given N in a simple closed network.
%
% J.Brooks
% 12 Oct 2012

addpath inst;
close all;
clear all;

% Economic Parameters
ratio = 2;
R = 10;
C = R/ratio;

% Routing probability matrix
P = [0 1; 1 0];

% Calculate visit ratio (solve traffic equations)
V = qnvisits(P);

% Build nodes
lf = 20;
mu = 1/2;
travel= qnmknode('-/g/inf', lf, 0);  % deterministic
server = qnmknode('m/m/m-fcfs', 1/mu);
%S = [2 3];                              % mean service time

% Analyze closed network
% U - utilization vector
% W - response time (W = Wq + S)
% Q - average number of customers (L = Lq + p)
% X - throughput

figure(1); title('Throughput');
figure(2); title('Net Economic Benefit Per Loop');
figure(3); title('Response Time');
figure(4); title('Number of Customers');
figure(5); title('Utilization');
for N = 1:2:40
    %    [U W Q X] = qnclosed(N, S, V);
    [U W Q X] = qnsolve("closed", N, {travel, server}, V );

    figure(1);
    plot(N, X(1), 'bs', 'markersize', 20), hold on;
    plot(N, X(2), 'r*', 'markersize', 20);
    lambda =min(roots([1, -(N+1)/lf - mu, N/lf*mu]));
    plot(N,lambda , 'k*', 'markersize', 10);
    
    figure(2);
    plot(N, R - C*sum(W), '*', 'markersize', 10), hold on;

    figure(3);
    plot(N, W(1), 'bs', 'markersize', 20), hold on;
    plot(N, W(2), 'r*', 'markersize', 20);
    plot(N, 1/(mu-lambda), 'k*', 'markersize', 15);

    figure(4);
    plot(N, Q(1), 'bs', 'markersize', 20), hold on;
    plot(N, Q(2), 'r*', 'markersize', 20);

    figure(5);
    plot(N, U(1), 'bs', 'markersize', 20), hold on;
    plot(N, U(2), 'r*', 'markersize', 20);

end

figure(1); title('Throughput');
legend({'Travel', 'Service', 'Approximation'}, 'location', 'northwest');
figure(2); title('Net Economic Benefit Per Loop');
legend({'Travel', 'Service'}, 'location', 'northwest');
figure(3); title('Response Time (Latency) ');
legend({'Travel', 'Service', 'Approximation'}, 'location', 'northwest');
figure(4); title('Number of Customers');
legend({'Travel', 'Service'}, 'location', 'northwest');
figure(5); title('Utilization');
legend({'Travel', 'Service'}, 'location', 'northwest');