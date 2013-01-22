close all;
clear all;

type = [1, repmat(2,1,6)];
mu = [20, repmat(5,1,6)];
visits = [6, 3, 2, 1, 3, 3, 3, 3, 4, 5, 1];
METHOD = 1;         % M/M/1

d = zeros(1,6);

site1 = [1, 5, 6, 7, 9, 10];
d(1) = sum(visits(site1));

site2 = [1, 5, 6, 10];
d(2) = sum(visits(site2));

site3 = [1, 5, 11];
d(3) = sum(visits(site3));

site4 = [1, 2, 3, 4, 8, 9, 10];
d(4) = sum(visits(site4));

site5 = [1, 2, 3, 7, 8, 9, 10];
d(5) = sum(visits(site5));

site6 = [1, 2, 6, 7, 8, 9, 10];
d(6) = sum(visits(site6));

% d = d*0.05;
N_RANGE = [2:20]
prob = zeros(length(N_RANGE), 6);
prob2 = zeros(length(N_RANGE), 6);
throughput = zeros(1, length(N_RANGE));
for scale = [0.001:0.002:0.1]
for N = N_RANGE
    d2 = d*scale*N/10;
    [u,w,q,x,flag] = optimalAssignment(mu,type,d2,N,METHOD);
    [Nvec, p] = routing(type, mu, d2, x, N, METHOD);
    
    prob(N-1,:) = p(1,2:end);
    prob2(N-1,:) = Nvec/N;
    throughput(N-1) = sum(x);
end

figure(1), plot3(repmat(scale,1,length(N_RANGE)),N_RANGE, prob), hold on;
title('Optimal Routing Probabilities');
zlabel('Percent Routed from Staging to Region');
ylabel('Number of Vehicles');
xlabel('Travel Scale');
legend({'1', '2', '3', '4', '5', '6'});
figure(2), plot3(repmat(scale,1,length(N_RANGE)),N_RANGE,throughput), hold on;
title('Expected Throughput Rate');
zlabel('Cubic Yards per Hour');
ylabel('Number of Vehicles');
xlabel('Travel Scale');
end