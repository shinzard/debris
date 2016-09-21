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
dispatchParams.Update = true            % update prices?

compParams.deltaA = -0.3;               % composition effect
compParams.deltaN = -2.2;               % composition effect

% DoE
reps = 50;                              % number of replications
duration = 60;                          % length of simulation (hours)

% Network topology
mus = [10.26, 6.62, 6, 6, 6, 6, 6, 6];  % for parallel servers,
                                        % this is 'base mu'
type = [1 1 2 2 2 2 2 2];               % NOTE: simulateQ assumes
                                        % that the central servers
                                        % are listed first (type=1)
% Travel times (hrs)
dij = [ 0.17 0.62; ...                  % from prices.xlsx       
        0.23 0.61; ...
        0.24 3.35; ...
        0.22 0.76; ...
        0.47 0.49; ...
        0.50 1.11 ];

% Remote --> Central routing probabilities (Debris type)
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

% Generate all prices from simplex
p1 = linspace(0,1,11);
p2 = linspace(0,1,11);
p3 = linspace(0,1,11);
p4 = linspace(0,1,11);
p5 = linspace(0,1,11);
p6 = linspace(0,1,11);
[p1, p2, p3, p4, p5, p6] = ndgrid(p1, p2, p3, p4, p5, p6);
bad = (p1+p2+p3+p4+p5+p6 ~= 1 | p1 == 0 | p2 == 0 | p3 == 0 | p4 == ...
       0 | p5 == 0 | p6 == 0);
p1(bad) = NaN;
p2(bad) = NaN;
p3(bad) = NaN;
p4(bad) = NaN;
p5(bad) = NaN;
p6(bad) = NaN;

keep = find(~isnan(p1));
numPrices = length(keep);

priceOpt = [p1(keep), p2(keep), p3(keep), p4(keep), p5(keep), p6(keep)];

% Starting allocation (uniform)
numTrucks = 24;
numDrts = sum(type==2);
N = repmat(numTrucks/numDrts, 1, numDrts);

run = 1;
throughput = zeros(1,numPrices*reps + 2*reps);
revenue = zeros(1,numPrices*reps + 2*reps);
inequity = zeros(1,numPrices*reps + 2*reps);
inequity2 = zeros(1,numPrices*reps + 2*reps);
inequity3 = zeros(1,numPrices*reps + 2*reps);
inequity4 = zeros(1,numPrices*reps + 2*reps);
T = zeros(1,numPrices*reps + reps);
priceStrategy = zeros(1,numPrices*reps + reps);

for temp = [0.73, 0.2]
    dispatchParams.T = temp            % level of rationality

    for p = 1:numPrices
        prices = priceOpt(p,:)         % prices
    
        for i = 1:reps
            fprintf('.');
            [u,w,q,x] =  simulateQ(N, mus, type, -1, dij, dispatchParams, ...
                                   compParams, pr, 60, 0, prices);
            
            throughput(run) = sum(x(find(type==1)));
            revenue(run) = sum(prices.*x(type==2));
            inequity(run) =  sqrt(sum(((x(type==2)./mean(x(type==2))-1).^2))/(length(x)*length(x)-1));% Flow
            inequity2(run) = sqrt(sum(((w(type==2)./mean(w(type==2))-1).^2))/(length(x)*length(x)-1));% Wait
            inequity3(run) = sqrt(sum(((q(type==2)./mean(q(type==2))-1).^2))/(length(x)*length(x)-1));% Length
            inequity4(run) = sqrt(sum(((u(type==2)./mean(u(type==2))-1).^2))/(length(x)*length(x)-1));% Utilization
            T(run) = temp;
            priceStrategy(run) = p;
            run = run + 1;
        end
    end
end

% Just run optimal prices...
prices = [1, 0.9345, 0.9515, 1.0875, 0.9438, 1.1404];% optimal
prices = prices/sum(prices);

for temp = [0.73, 0.2]
    dispatchParams.T = temp            % level of rationality
    for i = 1:reps
        fprintf('.');
        [u,w,q,x] =  simulateQ(N, mus, type, -1, dij, dispatchParams, ...
                               compParams, pr, 60, 0, prices);
        
        throughput(run) = sum(x(find(type==1)));
        revenue(run) = sum(prices.*x(type==2));
        inequity(run) =  sqrt(sum(((x(type==2)./mean(x(type==2))-1).^2))/(length(x)*length(x)-1));% Flow
        inequity2(run) = sqrt(sum(((w(type==2)./mean(w(type==2))-1).^2))/(length(x)*length(x)-1));% Wait
        inequity3(run) = sqrt(sum(((q(type==2)./mean(q(type==2))-1).^2))/(length(x)*length(x)-1));% Length
        inequity4(run) = sqrt(sum(((u(type==2)./mean(u(type==2))-1).^2))/(length(x)*length(x)-1));% Utilization
        T(run) = temp;
        priceStrategy(run) = numPrices+1;
        run = run + 1;
    end
end

% Now just fixed optimal partition...
dispatchParams.method = 'Optimal'      % {'Fixed', 'Bounded', 'Optimal'}

for i = 1:reps
    fprintf('.');
    [u,w,q,x] =  simulateQ(N, mus, type, -1, dij, dispatchParams, ...
                           compParams, pr, 60, 0, prices);
    
    throughput(run) = sum(x(find(type==1)));
    revenue(run) = sum(prices.*x(type==2));
    inequity(run) =  sqrt(sum(((x(type==2)./mean(x(type==2))-1).^2))/(length(x)*length(x)-1));% Flow
    inequity2(run) = sqrt(sum(((w(type==2)./mean(w(type==2))-1).^2))/(length(x)*length(x)-1));% Wait
    inequity3(run) = sqrt(sum(((q(type==2)./mean(q(type==2))-1).^2))/(length(x)*length(x)-1));% Length
    inequity4(run) = sqrt(sum(((u(type==2)./mean(u(type==2))-1).^2))/(length(x)*length(x)-1));% Utilization
    T(run) = NaN;
    priceStrategy(run) = numPrices+2;
    run = run + 1;
end

% PRINT!
% h = zeros(2,numPrices);
% for i = 1:numPrices
%     h(1,i) = mean(throughput(priceStrategy == i));
%     h(2,i) = 1.96*std(throughput(priceStrategy == i))/sqrt(reps);
% end

% figure, errorbar(h(1,:)/h(1,1), h(2,:)/h(1,1), 'b', 'LineWidth', ...
%                  0.75);

% h1 = mean(throughput(T==0.73 & priceStrategy==1));
% h2 = mean(throughput(T==0.73 & priceStrategy==2));
% h3 = mean(throughput(T==0.73 & priceStrategy==3));
% h4 = mean(throughput(T==0.73 & priceStrategy==4));
% l1 = mean(throughput(T==0.2 & priceStrategy==1));
% l2 = mean(throughput(T==0.2 & priceStrategy==2));
% l3 = mean(throughput(T==0.2 & priceStrategy==3));
% l4 = mean(throughput(T==0.2 & priceStrategy==4));

% opt = mean(throughput(priceStrategy == 5));

% h1_ci = 1.96*std(throughput(T==0.73 & priceStrategy==1))/sqrt(reps);
% h2_ci = 1.96*std(throughput(T==0.73 & priceStrategy==2))/sqrt(reps);
% h3_ci = 1.96*std(throughput(T==0.73 & priceStrategy==3))/sqrt(reps);
% h4_ci = 1.96*std(throughput(T==0.73 & priceStrategy==4))/sqrt(reps);
% l1_ci = 1.96*std(throughput(T==0.2 & priceStrategy==1))/sqrt(reps);
% l2_ci = 1.96*std(throughput(T==0.2 & priceStrategy==2))/sqrt(reps);
% l3_ci = 1.96*std(throughput(T==0.2 & priceStrategy==3))/sqrt(reps);
% l4_ci = 1.96*std(throughput(T==0.2 & priceStrategy==4))/sqrt(reps);

% opt_ci = 1.96*std(throughput(priceStrategy==5))/sqrt(reps);

% figure;
% errorbar([h1, h2, h3, h4, opt]/h1, [h1_ci, h2_ci, h3_ci, h4_ci, opt_ci]/h1, ...
%          'b', 'LineWidth', 0.75);
% hold on;
% errorbar([l1, l2, l3, l4]/h1, [l1_ci, l2_ci, l3_ci, l4_ci]/h1, ...
%          'b--', 'LineWidth', 0.75);

% set(gca,'FontSize', 15, 'XTick', [1, 2, 3, 4, 5], 'XTickLabel', ...
%         {'Uniform', 'OptimalPrice','Time', 'Dist', 'OptimalFixed'} , 'LineWidth',1)
% legend({'High T', 'Low T'}, 'FontSize', 12);
% %legend boxoff;

title('Mean Throughput');
ylabel('Percent of Uniform');
xlabel('Pricing Strategy');
%axis([0.5, 3.5, 0.7, 1.05])
