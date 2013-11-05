% This script is used to generate optimal assignments for the Fire Island
% case study (Hurricane Sandy)
%
% J.Brooks
%
% update 11 Feb 2013: reflects ARENA model travel times and routing
% (different only for team #3)
%
% update 20 Feb 2013: modified mu vector to have different pickup rates


close all;
clear all;

% Parameters
type = [1, repmat(2,1,6)];
METHOD = 1;         % M/M/1
N_RANGE = [6:6:60];
d = zeros(1,6);
mu1range = [5:5:40];
mu2range = [1:10];
scale = 10;

% Output matrices
prob = zeros(length(N_RANGE), 6);
prob2 = zeros(length(N_RANGE), 6);
throughput = zeros(1, length(N_RANGE));
numUniform = zeros(length(mu1range),length(mu2range));

% Set to either 'A' or 'B':
loc = 'A';

% TDSR Location A
if loc == 'A'
    visits = [0.039, 0.039, 0.017, 0.017, 0.011, 0.017, 0.017, ...
              0.039, 0.017, 0.011];
% TDSR Location B
elseif loc == 'B'
    visits = [0.006, NaN, 0.017, 0.017, 0.033, 0.017, 0.017, ...
              0.028, 0.017, 0.006];
end

% Block 1
if loc == 'A'
    site1 = [1, 5, 6, 7, 9, 10];
elseif loc == 'B'
    site1 = [1, 8];
end
d(1) = sum(visits(site1));

% Block 2
if loc == 'A'
    site2 = [1, 5, 6, 10];
elseif loc == 'B'
    site2 = [1, 7, 8, 9];
end
d(2) = sum(visits(site2));

% Block 3
if loc == 'A'
    site3 = [6, 7, 8, 9, 10];
elseif loc == 'B'
    site3 = [1, 6, 7, 8, 9, 10];
end
d(3) = sum(visits(site3));

% Block 4
if loc == 'A'
    site4 = [1, 2, 3, 4, 8, 9, 10];
elseif loc == 'B'
    site4 = [1, 5];
end
d(4) = sum(visits(site4));

% Block 5
if loc == 'A'
    site5 = [1, 2, 3, 7, 8, 9, 10];
elseif loc == 'B'
    site5 = [1, 4, 5, 7];
end
d(5) = sum(visits(site5));

% Block 6
if loc == 'A'
    site6 = [1, 2, 6, 7, 8, 9, 10];
elseif loc == 'B'
    site6 = [1, 3, 4, 5, 6, 7];
end
d(6) = sum(visits(site6));

% Particular case:
mu1range = 10
mu2range = 10

i = 1;
j = 1;
for mu1 = mu1range
    j = 1;
    for mu2 = mu2range
        %mu = [mu1, repmat(mu2,1,6)];
        
        % 20 Feb 2013: Latest different pickup-service rates
        mu = [1/0.05, 1/0.4, 1/0.3, 1/0.2, 1/0.6, 1/0.1, 1/0.5];
        idx = 1;
        for N = N_RANGE
            d2 = d*scale;
            [u,w,q,x,flag] = optimalAssignment(mu,type,d2,N,METHOD);
            [Nvec, p] = routing(type, mu, d2, x, N, METHOD);
    
            prob(idx,:) = p(1,2:end);
            prob2(idx,:) = Nvec/N;
            throughput(idx) = sum(x);
            idx = idx + 1;
        end

        numUniform(i,j) = sum(all(prob2 == 1/6,2));

        j = j + 1;
    end
    i = i + 1;
end

%figure(1), plot3(repmat(scale,1,length(N_RANGE)),N_RANGE, prob), hold on;
%title('Optimal Routing Probabilities');
%zlabel('Percent Routed from Staging to Region');
%ylabel('Number of Vehicles');
%xlabel('Travel Scale');
%legend({'1', '2', '3', '4', '5', '6'});
%figure(2), plot3(repmat(scale,1,length(N_RANGE)),N_RANGE,throughput), hold on;
%title('Expected Throughput Rate');
%zlabel('Cubic Yards per Hour');
%ylabel('Number of Vehicles');
%xlabel('Travel Scale');
%end

%figure, plot(N_RANGE, prob2);

for i = 1:10
prob2(i,:)*(6*i)
end