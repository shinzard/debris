function [Nvec, p] = routing(type, mus, d, f , N, METHOD)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    chains = find(type == 1);
    cycles = find(type == 2);
    numChains = length(chains);
    numCycles= length(cycles);
    decisionVars = length(f);
    numStations = length(mus);
    centralFlow = zeros(1,numChains);
    cycleFlow = zeros(1,numCycles);
    centralL = zeros(1,numChains); 
    cycleL = zeros(1,numCycles); 

    travelL = f.*d;

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
    
    p = zeros(numStations);
    for j = 1:numChains
        for k = 1:numCycles
            idx = (j-1)*numCycles + k;
            p(j, numChains+k) = f(idx)/centralFlow(j);
            p(numChains+k, j) = f(idx)/cycleFlow(k);
        end
    end

end

