% Uncomment for execution on morse
%addpath('/usr/share/octave/packages/3.2/queueing-1.1.1');
close all;
clear all;
more off;

% These need to match the generated data from optimalAssignment
load rand_2_6_20nov2012

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

V = qnvisits(P);

% Results vectors
u = zeros(NUM_TESTS,numStations);
w = zeros(NUM_TESTS,numStations);
q = zeros(NUM_TESTS,numStations);
x = zeros(NUM_TESTS,numStations);
n1 = zeros(NUM_TESTS,1);

test = 1;


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
    
    for j=1:numChains
        centralFlow(j) = sum(f((j-1)*numCycles + 1:j* ...
                               numCycles));
        centralL(j) = centralFlow(j)/(mus(chains(j))-centralFlow(j));
    end
    
    for j = 1:numCycles
        cycleFlow(j) = sum(f(j:numCycles:end));
        cycleL(j) = cycleFlow(j)/(mus(cycles(j))-cycleFlow(j));
    end

    % Rounding  (pure strategy)
    tmp = 1;
    for j = 1:numChains
        for k = 1:numCycles
            Nvec(tmp) = floor(travelL(tmp) + centralL(j)*(f(tmp)/ ...
                                                          centralFlow(j)) ...
                              + cycleL(k)*(f(tmp)/cycleFlow(k)) ); 
            tmp = tmp + 1;
        end
    end
    
    remaining = N - sum(Nvec);
    [x,idx] = sort(f - Nvec, 'descend');
    Nvec(idx(1:remaining)) = Nvec(idx(1:remaining)) + 1;
    
    if sum(Nvec)~= N
        Nvec
        sum(Nvec)
        remaining
        idx
        f
        N
        error('bad rounding!!');
    end
    
    % compare pure assignment
    %    [ut,wt,qt,xt] = qnclosed(Nvec,S,V);
    [ut,wt,qt,xt] = queueSim2(Nvec,S,P,500,100);

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