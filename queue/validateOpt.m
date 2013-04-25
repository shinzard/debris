% Uncomment for execution on morse
addpath('/usr/share/octave/packages/3.2/queueing-1.1.1');
close all;
clear all;
more off;

% These need to match the generated data from optimalAssignment
load data/rand_2_6_20nov2012_fewerN_mm1k
METHOD = 2;

%load data/rand_2_6_20nov2012_fewerN_finiteSource
%METHOD = 3;

%load data/rand_2_6_20nov2012
%METHOD = 1;

%type = [1 1 2 2 2 2 2 2]; IN MAT FILE
%NUM_TESTS = 100; IN MAT FILE

cycles = find(type == 2);
chains = find(type == 1);

numCycles = length(cycles);
numChains = length(chains);
numStations = length(type);
decisionVars = numCycles*numChains;

throughputSIM = zeros(1,NUM_TESTS*decisionVars);
waitTimesSIM = zeros(1,NUM_TESTS*numStations);
utilizationsSIM = zeros(1,NUM_TESTS*numStations);
lengthsSIM = zeros(1,NUM_TESTS*numStations);

% System configuration (pure assignment)
P = zeros(decisionVars,numStations,decisionVars,numStations);
tmp = 1;
for i = 1:numChains
    for j = 1:numCycles
        P(tmp, i, tmp ,numChains+j) = 1;
        P(tmp, numChains+j, tmp, i) = 1;
        tmp = tmp + 1;
    end
end

% Results vectors
u = zeros(NUM_TESTS,numStations);
w = zeros(NUM_TESTS,numStations);
q = zeros(NUM_TESTS,numStations);
x = zeros(NUM_TESTS,numStations);
n1 = zeros(NUM_TESTS,1);

disp(sprintf('Num \t Max \t Opt \t Sim'));

for i = 1:NUM_TESTS
    % Extract continuous assignment
    f = flows((i-1)*decisionVars + 1 : i*decisionVars);
    d = travelDelay((i-1)*decisionVars + 1 : i*decisionVars);
    mus = serviceRates((i-1)*numStations + 1 : i*numStations);
    N = entities(i);
    
    S = repmat(1./mus,decisionVars,1);

    centralFlow = zeros(1,numChains);
    cycleFlow = zeros(1,numCycles);
    centralL = zeros(1,numChains); 
    cycleL = zeros(1,numCycles); 

    travelL = f.*d;
    
    %    for j=1:numChains
    %        centralFlow(j) = sum(f((j-1)*numCycles + 1:j* ...
    %                               numCycles));
    %        centralL(j) = centralFlow(j)/(mus(chains(j))-centralFlow(j));
    %    end
    
    %    for j = 1:numCycles
    %        cycleFlow(j) = sum(f(j:numCycles:end));
    %        cycleL(j) = cycleFlow(j)/(mus(cycles(j))-cycleFlow(j));
    %    end

    for j=1:numChains
        centralFlow(j) = sum(f((j-1)*numCycles + 1:j* ...
                               numCycles));
        switch(METHOD)
          case 1,                       % M/M/1
            delta = centralFlow(j)/(mus(chains(j))-centralFlow(j));
            pk = 0;
            
          case 2,                       % M/M/1/N
            rho = centralFlow(j)/mus(chains(j));
            if rho < 1
                pk = ( rho^N*(1-rho) )/( 1 - rho^(N+1) );
                Lq = rho/(1-rho) - rho*(N*rho^N + 1)/(1-rho^(N+1));
            else
                pk = 1/(N+1);
                Lq = N*(N+1)/( 2*(N+1) );
            end
            delta = Lq + rho*(1-pk);
          case 3,                       % finite-source
            rho = centralFlow(j)/mus(chains(j));
            a = zeros(1,N);
            for k = 1:N
                a(k) = nchoosek(N,k)*factorial(k)*rho^k;
            end
            p0 = 1/(1+sum(a));
            delta = p0*sum([1:N].*a);
        end
        centralL(j) = delta; %min(max(0,delta), N);
    end
    
    for j = 1:numCycles
        cycleFlow(j) = sum(f(j:numCycles:end));

        switch(METHOD)
          case 1,                       % M/M/1
            delta = cycleFlow(j)/(mus(cycles(j))-cycleFlow(j));
            pk = 0;
            
          case 2,                       % M/M/1/N
            rho = cycleFlow(j)/mus(cycles(j));
            if rho < 1
                pk = ( rho^N*(1-rho) )/( 1 - rho^(N+1) );
                Lq = rho/(1-rho) - rho*(N*rho^N + 1)/(1-rho^(N+1));
            else
                pk = 1/(N+1);
                Lq = N*(N+1)/( 2*(N+1) );
            end
            delta = Lq + rho*(1-pk);
          case 3,                       % finite-source
            rho = cycleFlow(j)/mus(cycles(j));
            a = zeros(1,N);
            for k = 1:N
                a(k) = nchoosek(N,k)*factorial(k)*rho^k;
            end
            p0 = 1/(1+sum(a));
            delta = p0*sum([1:N].*a);
        end
        
        cycleL(j) = delta; %min(max(0,delta), N);
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

    NvecOrig = Nvec;

    Nvec = floor(Nvec);
    remaining = N - sum(Nvec);
    [x,idx] = sort(NvecOrig - Nvec, 'descend');
    Nvec(idx(1:remaining)) = Nvec(idx(1:remaining)) + 1;
    
    if sum(Nvec)~= N
        Nvec
        sum(Nvec)
        remaining
        idx
        f
        N
        sum(centralL) + sum(cycleL)
        error('bad rounding!!');
    end
    
    % pure assignment
    %    [ut,wt,qt,xt] = queueSim2(Nvec,S,P,500,100);

    %   Outputs:
    %     U - utilization vector
    %     W - response time (W = Wq + S)
    %     Q - average number of customers (L = Lq + p)
    %     X - throughput
    %
    %   Inputs:
    %     N - vector of entities by routing class (length = # parallel cycles)
    %     S - vector of mean service times
    %     P - routing probability vector (0 if determinstic routing by
    %     class, -1 for rational routing)


    % Pure assignment #2
    %    [ut2,wt2,qt2,xt2] = queueSim2(Nvec,S,P,500,100);

    % Probabilistic Assignment
    P2 = zeros(numStations);
    for j = 1:numChains
        for k = 1:numCycles
            idx = (j-1)*numCycles + k;
            P2(j, numChains+k) = f(idx)/centralFlow(j);
            P2(numChains+k, j) = f(idx)/cycleFlow(k);
        end
    end

    [ut,wt,qt,xt] = queueSim2(N,1./mus,P2,500,100);

    % for finite source only!!! (not sure why doesn't work when use
    % cycles instad of centrals...)
    if METHOD == 3
        f = sum(centralFlow.*(N-centralL));
    end
    
    disp(sprintf('%d \t %2.2f \t %2.2f \t %2.2f', i, min(sum(mus(chains)), ...
                                                    sum(mus(cycles))), ...
                 sum(f), sum(sum(xt(:,chains)))));
    % stash raw data
    throughputSIM((i-1)*numStations + 1 : i*numStations) = xt;
    waitTimesSIM((i-1)*numStations + 1 : i*numStations) = wt;
    utilizationsSIM((i-1)*numStations + 1 : i*numStations) = ut;
    lengthsSIM((i-1)*numStations + 1 : i*numStations) = qt;

end

save validationData throughputSIM waitTimesSIM utilizationsSIM lengthsSIM;