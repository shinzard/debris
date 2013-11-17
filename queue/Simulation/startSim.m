% This is the calling script to run the simulation experiments. 
% 
% 4 Nov 2013
% J.Brooks
% 
% Last Modified: 4 Nov 2013
% 
clear all;

% Parameters
dispatchParams.method = 'Bounded';      % {'Fixed', 'Bounded'}
dispatchParams.T = 0.73;                % Boltzmann temperature parameter
dispatchParams.inter = 0.44;            % hours
dispatchParams.forgetFactor = 0.9;      % EWMA forgetting factor

compParams.deltaA = 1;                  % composition effect
compParams.deltaN = 1;                  % composition effect

% DoE
reps = 100;                             % number of replications
duration = 60;                          % hours

% Network topology
mus = [10.26, 6.62, 6, 6, 6, 6, 6, 6];
type = [1 1 2 2 2 2 2 2];               % NOTE: simulateQ assumes
                                        % that the central servers
                                        % are listed first (type=1)
% Travel times (hrs)
%dij = [ 0.17 0.62; ...                  % from prices.xlsx       
%        0.23 0.61; ...
%        0.24 3.35; ...
%        0.22 0.76; ...
%        0.47 0.49; ...
%        0.50 1.11 ];
dij = [ 0.50 1.11; ...                  % from ARENA Model
                                        % (reversed row order)
        0.47 0.49; ...
        0.22 0.76; ...
        0.24 3.35; ...        
        0.23 0.61; ...
        0.17 0.62];



% Remote --> Central routing probabilities (Debris type)
pr = [ 0.33 0.67; ...
       0.60 0.40; ...
       0.92 0.08; ...
       0.20 0.80; ...
       0.92 0.08; ...
       0.43 0.57 ];

% Pricing strategies
priceOpt = [1 1 1 1 1 1; ...              % uniform
            0.85, 0.69, 0.91, 1.18, 0.85, 1.53; ...% time
            1.70, 0.74, 0.29, 1.93, 0.28, 1.06];% distance

% Starting allocation
numTrucks = 24;
numDrts = sum(type==2);
N = repmat(numTrucks/numDrts, 1, numDrts);

run = 1;

for temp = [0.73, 0.2]
    dispatchParams.T = temp            % level of rationality

    for p = 1:3
        prices = priceOpt(p,:)         % prices
    
        for i = 1:reps
            fprintf('.');
            [u,w,q,x] =  simulateQ(N, mus, type, -1, dij, dispatchParams, ...
                                   compParams, pr, 60, 0, prices);
            
            throughput(run) = sum(x(find(type==1)));
            T(run) = temp;
            priceStrategy(run) = p;
            run = run + 1;
        end
    end
end


h1 = mean(throughput(T==0.73 & priceStrategy==1));
h2 = mean(throughput(T==0.73 & priceStrategy==2));
h3 = mean(throughput(T==0.73 & priceStrategy==3));
l1 = mean(throughput(T==0.2 & priceStrategy==1));
l2 = mean(throughput(T==0.2 & priceStrategy==2));
l3 = mean(throughput(T==0.2 & priceStrategy==3));

h1_ci = 1.96*std(throughput(T==0.73 & priceStrategy==1))/sqrt(reps);
h2_ci = 1.96*std(throughput(T==0.73 & priceStrategy==2))/sqrt(reps);
h3_ci = 1.96*std(throughput(T==0.73 & priceStrategy==3))/sqrt(reps);
l1_ci = 1.96*std(throughput(T==0.2 & priceStrategy==1))/sqrt(reps);
l2_ci = 1.96*std(throughput(T==0.2 & priceStrategy==2))/sqrt(reps);
l3_ci = 1.96*std(throughput(T==0.2 & priceStrategy==3))/sqrt(reps);

figure, plot([l1, l2, l3]/h1, 'b--');
hold on;
plot([h1, h2, h3]/h1, 'b');

figure;
errorbar([h1, h2, h3]/h1, [h1_ci, h2_ci, h3_ci]/h1, ...
         'b', 'LineWidth', 0.75);
hold on;
errorbar([l1, l2, l3]/h1, [l1_ci, l2_ci, l3_ci]/h1, ...
         'b--', 'LineWidth', 0.75);

set(gca,'FontSize', 15, 'XTick', [1, 2, 3], 'XTickLabel', ...
        {'Uniform', 'Time', 'Dist'} , 'LineWidth',1)
legend({'High T', 'Low T'}, 'FontSize', 12);
%legend boxoff;
title('Mean Throughput');
ylabel('Percent of Uniform');
xlabel('Pricing Strategy');
axis([0.5, 3.5, 0.8, 1.05])
