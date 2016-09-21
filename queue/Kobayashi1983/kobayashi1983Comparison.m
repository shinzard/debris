clear all;

addpath('../');
% Number of simulations to run
SAMPLES = 100;

mus = [4 2 1 0.5];
type = [1 2 2 2];

N = 8;
METHOD = 4;

[u,w,q,x,e]=optimalAssignment(mus, type, 0, N, METHOD)
[Nvec, p] = routing(type, mus, 0, x, N, METHOD)

P = [0,x./sum(x); 1 0 0 0; 1 0 0 0;1 0 0 0];

% from Table V on pg. 305
% for N = 5
%P2 = [0,0.693, 0.239, 0.068; 1 0 0 0; 1 0 0 0;1 0 0 0];

% For N = 8
P2 = [0,0.644, 0.260, 0.097; 1 0 0 0; 1 0 0 0;1 0 0 0];

kResults = zeros(1,SAMPLES);
results = zeros(1,SAMPLES);
for i = 1:SAMPLES
    disp(sprintf('%d',i));
    [u,w,q,x]=queueSim2(N,1./mus,P,5000,1000);
    results(i) = x(1);
    
    %    [u,w,q,x]=queueSim2(N,1./mus,P2,5000,1000);
    %    kResults(i) = x(1);

    [u,w,q,x]=queueSim2(Nvec,1./mus,0,5000,1000);
    FixedResults(i) = x(1);
    
    %disp(sprintf('%d: %2.2f, %2.2f, %2.2f',i,results(i), kResults(i), FixedResults(i)));
end

beep;

figure, ax(1)=subplot(311), hist(results), title('My Results');
ax(2)=subplot(312), hist(FixedResults), title('Fixed');
ax(3)=subplot(313), hist(kResults), title('Other Results');

linkaxes(ax, 'x');
subplot(311), axis([2.3, 2.95, 0, 25]);
title('Histogram of FPM Average Throughput');
subplot(312), axis([2.3, 2.95, 0, 25]);
title('Histogram of Customer-Deterministic Average Throughput');
subplot(313), axis([2.3, 2.95, 0, 25]);
title('Histogram of MVA Average Throughput');
xlabel('Throughput');


mean(results)
mean(FixedResults)
mean(kResults)

%save kobayashiCompare_withFixed_9June2013_mm1