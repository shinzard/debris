% This script is intended to determine the approximate delay
% function (and flow) for a given N in a double closed network.
%
% J.Brooks
% 17 Oct 2012

addpath inst;
close all;
clear all;

% Economic Parameters
ratio = 25;
R = 1;
C = R/ratio;

% Routing probability matrix
P = [0 0.5 0.5; 1 0 0; 1 0 0];
P1 = [0 0.8 0.2; 1 0 0; 1 0 0];
P2 = [0 0.2 0.8; 1 0 0; 1 0 0];

% Calculate visit ratio (solve traffic equations)
V = qnvisits(P);
V1 = qnvisits(P1);
V2 = qnvisits(P2);

% Build nodes
mu2 = 2;
mu = 1;
mu3 = 1.5;
%travel= qnmknode('-/g/inf', lf, 0);  % deterministic
Q1s = qnmknode('m/m/m-fcfs', 1/mu);
Q2s = qnmknode('m/m/m-fcfs', 1/mu);
TDSR = qnmknode('m/m/m-fcfs', 1/mu2);
%TDSR= qnmknode('-/g/inf', 1/mu2, 0);  % deterministic

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
    [U W Q X] = qnsolve("closed", N, {TDSR, Q1s, Q2s}, V);
    [U1 W1 Q1 X1] = qnsolve("closed", N, {TDSR, Q1s, Q2s}, V1);
    [U2 W2 Q2 X2] = qnsolve("closed", N, {TDSR, Q1s, Q2s}, V2);

    figure(1);
    plot(N, X(1), 'bs', 'markersize', 20), hold on;
    plot(N, X(2), 'r*', 'markersize', 20);
    %    lambda =min(roots([1, -(N+1)/lf - mu, N/lf*mu]));
    %    plot(N,lambda , 'k*', 'markersize', 10);
    
    figure(2);
    plot(N, R - C*(W(1)+0.5*(W(2)+W(3))), 'r*', 'markersize', 10);
    hold on;
    plot(N, R - C*(W1(1)+0.5*(W1(2)+W1(3))), 'b*', 'markersize', 10);
    plot(N, R - C*(W2(1)+0.5*(W2(2)+W2(3))), 'k*', 'markersize', 10);
    
    plot(N, N*(R - C*(W(1)+0.5*(W(2)+W(3)))), 'rs', 'markersize', ...
         10)
    plot(N, N*(R - C*(W1(1)+0.5*(W1(2)+W1(3)))), 'bs', 'markersize', 10)
    plot(N, N*(R - C*(W2(1)+0.5*(W2(2)+W2(3)))), 'ks', 'markersize', 10)

    figure(3);
    plot(N, W(1), 'bs', 'markersize', 20), hold on;
    plot(N, W(2), 'r*', 'markersize', 20);
    %    plot(N, 1/(mu-lambda), 'k*', 'markersize', 15);

    figure(4);
    plot(N, Q(1), 'bs', 'markersize', 20), hold on;
    plot(N, Q(2), 'r*', 'markersize', 20);

    figure(5);
    plot(N, U(1), 'bs', 'markersize', 20), hold on;
    plot(N, U(2), 'r*', 'markersize', 20);

end

figure(1); title('Throughput');
legend({'Central', 'Service'}, 'location', 'northwest');
figure(2); title('Net Economic Benefit Per Loop');
legend({'Individual', 'Social'}, 'location', 'northwest');
figure(3); title('Response Time (Latency) ');
legend({'Central', 'Service'}, 'location', 'northwest');
figure(4); title('Number of Customers');
legend({'Central', 'Service'}, 'location', 'northwest');
figure(5); title('Utilization');
legend({'Central', 'Service'}, 'location', 'northwest');