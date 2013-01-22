
% Number of simulations to run
SAMPLES = 100;

mus = [4 2 1 0.5];
type = [1 2 2 2];


[u,w,q,x,e]=optimalAssignment(mus,type,0,5,2);

P = [0,x./sum(x); 1 0 0 0; 1 0 0 0;1 0 0 0];

% from Table V on pg. 305
P2 = [0,0.693, 0.239, 0.068; 1 0 0 0; 1 0 0 0;1 0 0 0];

kResults = zeros(1,SAMPLES);
results = zeros(1,SAMPLES);
for i = 1:SAMPLES
    [u,w,q,x]=queueSim2(5,1./mus,P,5000,1000);
    results(i) = x(1);
    
    [u,w,q,x]=queueSim2(5,1./mus,P2,5000,1000);
    kResults(i) = x(1);
end

figure, subplot(211), hist(results), title('My Results');
subplot(212), hist(kResults), title('Other Results');