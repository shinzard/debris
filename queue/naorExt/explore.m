function [tollRate,ns, indivRate, thru, socialRate, socialThru] = explore(C,R,T,mu, PLOT)

if nargin < 5
    PLOT = 0;
end

more off;

addpath '~/Documents/MATLAB/queueing/inst'
%close all;
%clear all;

%% ORIGINAL
%R = 10;
%C = 10;

%mu = [10.26, 6.62, 6, 6, 6, 6, 6, 6];    % system from WSC paper...
%mu = [10.26, 6.62, 6, 7, 8, 9, 10, 11];
%P = zeros(8,8);
%P(1,3) = P(1,4) = P(1,5) = P(1,6) = P(1,7) = P(1,8) = 1/6;
%P(2,3) = P(2,4) = P(2,5) = P(2,6) = P(2,7) = P(2,8) = 1/6;
%P(3,1) = 0.33; P(3,2) = 0.67;
%P(4,1) = 0.60; P(4,2) = 0.40;
%P(5,1) = 0.92; P(5,2) = 0.08;
%P(6,1) = 0.20; P(6,2) = 0.80;
%P(7,1) = 0.92; P(7,2) = 0.08;
%P(8,1) = 0.43; P(8,2) = 0.57;

%V = qnvisits(P);

%total = 20;

%i = 1;
%U = zeros(8,total);
%W = zeros(8,total);
%Q = zeros(8,total);
%X = zeros(8,total);

%for N = 1:total
%    disp(sprintf('N: %d', N));
%    [U(:,i), W(:,i), Q(:,i), X(:,i)] = qnclosed(N,1./mu,V);
%     i = i + 1;
% end

% throughput = sum(X(1:2,:));

% wait = sum(P)./[6 6 2 2 2 2 2 2]*W

% socialProfit = (R - wait*C).*[1:total]

% plot([1:total], wait*C)
% hold on;
% plot([1:total], socialProfit, 'g');
% plot([1 total], [R R], 'r--');

%% SIMPLE SYSTEM FOR PAPER...
%...TODO: make consistent?? kind of a lot of work for
% this paper...
%C = 1;
%R = 5;
%T = 2;
%p = 0.1;                                % open feedback probability

%mu = [4 6];
Prob = [0 1; 1 0];                         % equivalent closed network

V = qnvisits(Prob);

%total = 150;
total = 25;

i = 1;
U = zeros(2,total);
W = zeros(2,total);
Q = zeros(2,total);
X = zeros(2,total);

for N = 1:total
    [U(:,i), W(:,i), Q(:,i), X(:,i)] = qnclosed(N,1./mu,V);
    i = i + 1;
    % TODO: add break when indiv loop profit negative...
end

throughput = X(1,:);
%netThroughput = throughput*(1-p);

wait = sum(W);

%indivProfit = (1/(1-p))*(R - wait*C) - T
indivProfit = (R - wait*C) - T;
%socialProfit = netThroughput.*(indivProfit + T) % include tolls
                                                 % here?? to be
                                                 % consistent with
                                                 % Naor, YES
socialProfit = throughput.*(indivProfit + T);

toll = throughput*T;

if PLOT
    figure;
    %    plot([1:total], 1/(1-p)*wait*C+T)
    plot([1:total], wait*C+T)
    hold on;
    set(gca, 'fontsize', 20);
    plot([1:total], socialProfit, 'g');
    %    plot([1 total], 1/(1-p).*[R R], 'r--');
    plot([1 total], [R R], 'r--');
    %plot([1:total], toll, 'k')
    %    plot([1:total], netThroughput, 'm')
    plot([1:total], throughput, 'm')
    legend({'Individual Cost', 'Social Benefit Rate', 'Indiv. Benefit', ...    %'Toll Rate', 
            'NetThru'}, 'location', 'northeast');
    title('Performance Measures');
    xlabel('Number of Customers');
    disp(sprintf('\t\t Num \t NetThru \t Social \t Indiv. \t `Tolls'));

    [val, idx] = max(socialProfit);
    disp(sprintf('Social Optimum: %d \t %f \t %f \t %f \t %f', idx, throughput(idx), ...
                 socialProfit(idx), indivProfit(idx), toll(idx)))

    idx = find(indivProfit>0, 1, 'last');
    disp(sprintf('Indiv. Optimum: %d \t %f \t %f \t %f \t %f', idx, throughput(idx), ...
                 socialProfit(idx), indivProfit(idx), toll(idx)));
end

% OUTPUTS
idx = find(indivProfit>0, 1, 'last');
if ~isempty(idx)
    tollRate = toll(idx);
    ns = idx;
    indivRate = indivProfit(idx)*throughput(idx);
    thru = throughput(idx);

    [val, idx] = max(socialProfit);
    socialRate = indivProfit(idx)*throughput(idx);
    socialThru = throughput(idx);
else
    tollRate = 0;
    ns = 0;
    indivRate = 0;
    thru = 0;
    socialRate = 0;
    socialThru = 0;
end
    
