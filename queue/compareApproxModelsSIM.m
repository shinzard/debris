load data/rand_2_6_20nov2012
load data/rand_2_6_20nov2012_mm1k
load data/rand_2_6_20nov2012_finiteSource

numStations = length(type);
numCycles = length(cycles);
numChains = length(chains);
numStations = length(type);
decisionVars = numCycles*numChains;

for i = 1:NUM_TESTS
    idx = [(i-1)*numStations + 1 : i*numStations];

    for j=1:numChains
        centralFlow(j) = sum(f((j-1)*numCycles + 1:j* ...
                               numCycles));
    end
    for j = 1:numCycles
        cycleFlow(j) = sum(f(j:numCycles:end));
    end

    % Rounding  (pure strategy)
    tmp = 1;
    Nvec = zeros(1,decisionVars);
    for j = 1:numChains
        for k = 1:numCycles
            if cycleFlow(k) > 0 && centralFlow(j) > 0
                Nvec(tmp) = travelL(tmp) + centralL(j)*(f(tmp)/ ...
                                                        centralFlow(j)) ...
                    + cycleL(k)*(f(tmp)/cycleFlow(k)); 
            end
            tmp = tmp + 1;
        end
    end

    % Probabilistic Assignment
    P2 = zeros(numStations);
    for j = 1:numChains
        for k = 1:numCycles
            idx = (j-1)*numCycles + k;
            P2(j, numChains+k) = f(idx)/centralFlow(j);
            P2(numChains+k, j) = f(idx)/cycleFlow(k);
        end
    end

    NvecOrig = Nvec;

    Nvec = floor(Nvec);
    remaining = N - sum(Nvec);
    [x,idx] = sort(NvecOrig - Nvec, 'descend');
    Nvec(idx(1:remaining)) = Nvec(idx(1:remaining)) + 1;


    %throughputSIM = zeros(1,NUM_TESTS*decisionVars);
    %waitTimesSIM = zeros(1,NUM_TESTS*numStations);
    %utilizationsSIM = zeros(1,NUM_TESTS*numStations);
    %lengthsSIM = zeros(1,NUM_TESTS*numStations);


% Finite-Source
[u,w,q,x,e]=optimalAssignment(mu,type,0,20,3)
figure(1), 
subplot(411), plot(u,'r.'), hold on;
subplot(412), plot(w,'r.'), hold on;
subplot(413), plot(q,'r.'), hold on;
subplot(414), plot(x,'r.'), hold on;

subplot(411), title('Utilization');
subplot(412), title('Wait Times');
subplot(413), title('Length')
subplot(414), title('Flows');

% M/M/1/N
[u,w,q,x,e]=optimalAssignment(mu,type,0,20,2);
subplot(411), plot(u,'g.');
subplot(412), plot(w,'g.');
subplot(413), plot(q,'g.');
subplot(414), plot(x,'g.');

% M/M/1
[u,w,q,x,e]=optimalAssignment(mu,type,0,20,1);
subplot(411), plot(u,'b.');
subplot(412), plot(w,'b.');
subplot(413), plot(q,'b.');
subplot(414), plot(x,'b.');

legend({'finite-source', 'mm1k', 'mm1'});