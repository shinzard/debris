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
dispatchParams.forgetFactor = 0.9;            % EWMA forgetting factor

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
dij = [0.47 0.47; ...
       0.38 0.38; ...
       0.50 0.50; ...
       0.65 0.65; ...
       0.47 0.47; ...
       0.85 0.85];

numDrts = sum(type==2);
numTrucks = 24;

% Initial distribution
N = repmat(numTrucks/numDrts, 1, numDrts);

% Remote --> Central routing probabilities (Debris type)
pr = [0.33 0.67; ...
      0.60 0.40; ...
      0.92 0.08; ...
      0.20 0.80; ...
      0.92 0.08; ...
      0.43 0.57];

% Starting allocation
assign = [4 4 4 4 4 4];

priceOpt = [1 1 1 1 1 1; ...              % uniform
            0.85, 0.69, 0.91, 1.18, 0.85, 1.53; ...% time
            1.7, 0.74, 0.29, 1.93, 0.28, 1.06];% distance

run = 1;
for temp = [0.73, 0.2]
    dispatchParams.T = temp            % level of rationality

    for p = 1:3
        prices = priceOpt(p,:)         % prices
    
        for i = 1:reps
            fprintf('.');
            [u,w,q,x] =  simulateQ(N, mus, type, -1, dij, dispatchParams, ...
                                   compParams, pr, 60, 0);
            
            throughput(run) = sum(x(find(type==1)));
            T(run) = temp;
            priceStrategy(run) = p;
            run = run + 1;
        end
    end
end