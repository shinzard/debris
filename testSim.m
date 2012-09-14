% This script simulates a simple debris situation
% 
% 19 October 2011
% J.Brooks

close all;
clear all;

% Parameters
SIM_TIME = 10*60;                       % Total simulation time in minutes
NUM_TRUCKS = 100;
NUM_SITES = 15;                          % number of loading sites (<100)
NUM_TDSR = 3;                           % number of TDSRs (<50)
MAX_RATE = 400;                         % maximum loader rate (cy/hr)
MIN_RATE = 200;                          % minimum loader rate
                                        % (cy/hr)
LOADER_VAR = 10;                        % loader rate variance
MIN_PERC_FULL = 0.40;                   % minimum percent full
MAX_CAP = 100;                           % maximum truck capacity
MIN_CAP = 20;                           % minimum truck capacity

% Initialize vehicles
loaderRate = randi([MIN_RATE, MAX_RATE], NUM_SITES, 1);
truckCapacity = randi([MIN_CAP, MAX_CAP], NUM_TRUCKS, 1);
truckLoad = zeros(1, NUM_TRUCKS);
truckStatus = ones(1, NUM_TRUCKS);
% States:
%   0    waiting at TDSR
%   1    waiting at loading site
%   2    enroute to TDSR
%   3    enroute to loading site

truckTimeStamp = zeros(1, NUM_TRUCKS);

% Initialize sites
siteAmount = zeros(1, NUM_SITES);

% Random initial assignments
truckTdsrAssignment = randi([1, NUM_TDSR], 1, NUM_TRUCKS);
truckLoadAssignment = randi([1, NUM_SITES], 1, NUM_TRUCKS);

% Storage data
amount = zeros(SIM_TIME, NUM_SITES);

% Run simulation
for t = 1:SIM_TIME
    
    % Update truck status as needed
    truckUpdate = find(truckTimeStamp == t);
    tdsr = find(truckStatus(truckUpdate) == 2);
    truckStatus(truckUpdate(tdsr)) = 0;
    load = find(truckStatus(truckUpdate) == 3);
    truckStatus(truckUpdate(load)) = 1;
    
    % Process loading sites
    for i = 1:NUM_SITES
        % find trucks waiting at site i
        trucksWaiting = find(truckStatus == 1 & truckLoadAssignment ...
                             == i);
        [tmp, idx] = min(truckTimeStamp(trucksWaiting));
        idx = trucksWaiting(idx);
        
        if trucksWaiting & (siteAmount(i) > truckCapacity(idx)*MIN_PERC_FULL)
            % Load truck
            amt = min(siteAmount(i), truckCapacity(idx));
            truckLoad(idx) = amt;

            truckStatus(idx) = 2;
            % Arrival time at TDSR
            truckTimeStamp(idx) = t + randi([40, 80], 1);
                        
            siteAmount(i) = siteAmount(i) - amt;
        else                            % continue sorting...
            siteAmount(i) = siteAmount(i) + loaderRate(i)/60 + LOADER_VAR/3600*randn(1);
        end
    end
    
    % Process TDSRs
    for j = 1:NUM_TDSR
        % find trucks waiting at tdsr j 
        trucksWaiting = find(truckStatus == 0 & truckTdsrAssignment ...
                             == j);
        [tmp, idx] = min(truckTimeStamp(trucksWaiting));
        idx = trucksWaiting(idx);
        
        if trucksWaiting
            truckLoad(idx) = 0;
            truckTimeStamp(idx) =  t + randi([40, 80], 1);
            truckStatus(idx) = 3;
        end
    end
    
    % Store time history
    amount(t, :) = siteAmount';
    
end

plot(amount)