function [u,w,q,x,exit] = optimalAssignmentWeights(mu,type,d,weights, ...
                                                   N,METHOD)

    chains = find(type == 1);
    numChains = length(chains);
    cycles = find(type == 2);
    numCycles = length(cycles);
    
    options = optimset('Algorithm', 'interior-point', 'TolX', 1e-12); %, ...
                       %'Display', 'iter');

    [x,fval, exit, output, lagrange] = ...
        fmincon(@(x)throughput(x,weights), repmat(0.0001,1,numChains*numCycles), [], ...
                [], [], [],repmat(0,1,numChains*numCycles),[], ...
                @(x)numEntities(x,mu,type,d,N,METHOD), options);
    
    centralFlow = zeros(1,numChains);
    cycleFlow = zeros(1,numCycles);
    
    centralL = zeros(1,numChains);
    cycleL = zeros(1,numCycles);
    travelL = sum(x.*d);
    
    centralW = zeros(1,numChains);
    cycleW = zeros(1,numCycles);
    
    %    for i=1:numChains
    %        centralFlow(i) = sum(x((i-1)*numCycles + 1:i*numCycles));
    %        centralL = centralL + 1/(mu(chains(i))-centralFlow(i));
    %    end
    
    %    for i = 1:numCycles
    %        cycleFlow(i) = sum(x(i:numCycles:end));
    %        cycleL = cycleL + 1/(mu(cycles(i))-cycleFlow(i));
    %    end
    
    for i=1:numChains
        % note that for finite-source, this is the *individual*
        % arrival rate, not the effective one...
        centralFlow(i) = sum(x((i-1)*numCycles + 1:i* ...
                               numCycles));
        switch(METHOD)
          case 1,                       % M/M/1
            delta = centralFlow(i)/(mu(chains(i))-centralFlow(i));
            pk = 0;
            
          case 2,                       % M/M/1/N
            rho = centralFlow(i)/mu(chains(i));
            if rho < 1
                pk = ( rho^N*(1-rho) )/( 1 - rho^(N+1) );
                Lq = rho/(1-rho) - rho*(N*rho^N + 1)/(1-rho^(N+1));
 
            else
                pk = 1/(N+1);
                Lq = N*(N+1)/( 2*(N+1) );
            end
            delta = Lq + rho*(1-pk);
            
          case 3,                     % finite source
            rho = centralFlow(i)/mu(chains(i));
            a = zeros(1,N);
            for j = 1:N
                a(j) = nchoosek(N,j)*factorial(j)*rho^j;
            end
            p0 = 1/(1+sum(a));
            pk = p0*a(end);
            delta = p0*sum([1:N].*a);
        end
        
        centralL(i) = delta; %min(max(0,delta), N);
        centralW(i) = centralL(i)/( centralFlow(i)*(1-pk) );
    end
    
    for i = 1:numCycles
        % note that for finite-source, this is the *individual*
        % arrival rate, not the effective one...
        cycleFlow(i) = sum(x(i:numCycles:end));

        switch(METHOD)
          case 1,                       % M/M/1
            delta = cycleFlow(i)/(mu(cycles(i))-cycleFlow(i));
            pk = 0;
            
          case 2,                       % M/M/1/N
            rho = cycleFlow(i)/mu(cycles(i));
            if rho < 1
                pk = ( rho^N*(1-rho) )/( 1 - rho^(N+1) );
                Lq = rho/(1-rho) - rho*(N*rho^N + 1)/(1-rho^(N+1));
 
            else
                pk = 1/(N+1);
                Lq = N*(N+1)/( 2*(N+1) );
            end
            delta = Lq + rho*(1-pk);

          case 3,                       % finite-source
            rho = cycleFlow(i)/mu(cycles(i));
            a = zeros(1,N);
            for j = 1:N
                a(j) = nchoosek(N,j)*factorial(j)*rho^j;
            end
            p0 = 1/(1+sum(a));
            pk = p0*a(end);
            delta = p0*sum([1:N].*a);
        end
        
        cycleL(i) = delta; %min(max(0,delta), N);
        cycleW(i) = cycleL(i)/( cycleFlow(i)*(1-pk) );
    end

    % Effective flows for finite source case
    if METHOD == 3
        centralFlow = centralFlow.*(N-centralL); 
        cycleFlow = cycleFlow.*(N-cycleL);
    end

    % Outputs
    u = [centralFlow./mu(chains), cycleFlow./mu(cycles)];
    switch(METHOD)
      case 1,
        w = [1./(mu(chains)-centralFlow), 1./(mu(cycles)- ...
                                              cycleFlow)];
      case 2,
        w = [centralW, cycleW];

      case 3,
        w = [centralW, cycleW];
    end
    q = [centralL, cycleL, travelL];
    
    disp(sprintf('System Throughput: %2.2f', sum(centralFlow)));

end


function T = throughput(lambda,weights)
    T = -sum(weights.*lambda);
end

function [C, Ceq] = numEntities(lambda, mu, type, d, N, METHOD)

    epsilon = 1e-6;                     % required to ensure
                                        % feasible solutions
    chains = find(type == 1);
    numChains = length(chains);
    cycles = find(type == 2);
    numCycles = length(cycles);
    
    centralFlow = zeros(1,numChains);
    cycleFlow = zeros(1,numCycles);
    
    centralL = 0;
    cycleL = 0;
    travelL = sum(lambda.*d);

    for i=1:numChains
        centralFlow(i) = sum(lambda((i-1)*numCycles + 1:i* ...
                                    numCycles));
        switch(METHOD)
          case 1,                       % M/M/1
            delta = centralFlow(i)/(mu(chains(i))-centralFlow(i));
            pk = 0;
            
          case 2,                       % M/M/1/N
            rho = centralFlow(i)/mu(chains(i));
            if rho < 1
                pk = ( rho^N*(1-rho) )/( 1 - rho^(N+1) );
                Lq = rho/(1-rho) - rho*(N*rho^N + 1)/(1-rho^(N+1));
 
            else
                pk = 1/(N+1);
                Lq = N*(N+1)/( 2*(N+1) );
            end
            delta = Lq + rho*(1-pk);

          case 3,                       % finite-source
            rho = centralFlow(i)/mu(chains(i));
            a = zeros(1,N);
            for j = 1:N
                a(j) = nchoosek(N,j)*factorial(j)*rho^j;
            end
            p0 = 1/(1+sum(a));
            delta = p0*sum([1:N].*a);
        end
        centralL = centralL + delta; %min(max(0,delta), N);
    end
    

    for i = 1:numCycles
        cycleFlow(i) = sum(lambda(i:numCycles:end));

        switch(METHOD)
          case 1,                       % M/M/1
            delta = cycleFlow(i)/(mu(cycles(i))-cycleFlow(i));
            pk = 0;
            
          case 2,                       % M/M/1/N
            rho = cycleFlow(i)/mu(cycles(i));
            if rho < 1
                pk = ( rho^N*(1-rho) )/( 1 - rho^(N+1) );
                Lq = rho/(1-rho) - rho*(N*rho^N + 1)/(1-rho^(N+1));
 
            else
                pk = 1/(N+1);
                Lq = N*(N+1)/( 2*(N+1) );
            end
            delta = Lq + rho*(1-pk);

          case 3,                       % finite-source
            rho = cycleFlow(i)/mu(cycles(i));
            a = zeros(1,N);
            for j = 1:N
                a(j) = nchoosek(N,j)*factorial(j)*rho^j;
            end
            p0 = 1/(1+sum(a));
            delta = p0*sum([1:N].*a);
        end
        
        cycleL = cycleL + delta; %min(max(0,delta), N);
    end
    
    C = [ centralL + cycleL + travelL - N , ...% total number
          centralFlow-mu(chains)+epsilon, ...  % stability of chains
          cycleFlow-mu(cycles)+epsilon];       % stability of cycles
    
    Ceq = [];
end



