function [u,w,q,x,lagrange] = optimalAssignment(mu,type,d,N)

    chains = find(type == 1);
    numChains = length(chains);
    cycles = find(type == 2);
    numCycles = length(cycles);
    
    options = optimset('Algorithm', 'interior-point');

    [x,fval, exit, output, lagrange] = ...
        fmincon(@(x)throughput(x), repmat(0,1,numChains*numCycles), [], ...
                [], [], [],repmat(0,1,numChains*numCycles),[], ...
                @(x)numEntities(x,mu,type,d,N), options);
    
        
    centralFlow = zeros(1,numChains);
    cycleFlow = zeros(1,numCycles);
    
    centralL = 0;
    cycleL = 0;
    travelL = 0;
    
    for i=1:numChains
        centralFlow(i) = sum(x((i-1)*numCycles + 1:i*numCycles));
        centralL = centralL + 1/(mu(chains(i))-centralFlow(i));
    end
    
    for i = 1:numCycles
        cycleFlow(i) = sum(x(i:numCycles:end));
        cycleL = cycleL + 1/(mu(cycles(i))-cycleFlow(i));
    end

    % Outputs
    u = [centralFlow./mu(chains), cycleFlow./mu(cycles)];
    w = [1./(mu(chains)-centralFlow), 1./(mu(cycles)-cycleFlow)];
    q = [centralFlow./(mu(chains)-centralFlow), cycleFlow./ ...
         (mu(cycles)-cycleFlow)];
    
    disp(sprintf('System Throughput: %2.2f', sum(centralFlow)));

end


function T = throughput(lambda)
    T = -sum(lambda);
end

function [C, Ceq] = numEntities(lambda, mu, type, d, N)
    chains = find(type == 1);
    numChains = length(chains);
    cycles = find(type == 2);
    numCycles = length(cycles);
    
    centralFlow = zeros(1,numChains);
    cycleFlow = zeros(1,numCycles);
    
    centralL = 0;
    cycleL = 0;
    travelL = 0;
    
    for i=1:numChains
        centralFlow(i) = sum(lambda((i-1)*numCycles + 1:i* ...
                                    numCycles));
        delta = centralFlow(i)/(mu(chains(i))-centralFlow(i));
        centralL = centralL + delta; %min(max(0,delta), N);
    end
    
    for i = 1:numCycles
        cycleFlow(i) = sum(lambda(i:numCycles:end));
        delta = cycleFlow(i)/(mu(cycles(i))-cycleFlow(i));
        cycleL = cycleL + delta; %min(max(0,delta), N);
    end
    
    C = [ centralL + cycleL + travelL - N , ... % total number
          centralFlow - mu(chains), ...         % stability of chains
          cycleFlow - mu(cycles) ];             % stability of cycles
    
    Ceq = [];
end



