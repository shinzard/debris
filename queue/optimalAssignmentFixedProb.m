% this allows for fixed probabilistic routing from parallel to
% central nodes via the input 'prob'
function [u,w,q,x,exit] = optimalAssignmentFixedProb(mu,type,d,weights,N,METHOD,remoteP)

    chains = find(type == 1);
    numChains = length(chains);
    cycles = find(type == 2);
    numCycles = length(cycles);
    
    options = optimset('Algorithm', 'interior-point', 'TolX', 1e-30, ...
                       'Display', 'off');

    [x,fval, exit, output, lagrange] = ...
        fmincon(@(x)throughput(x,weights), repmat(0.0001,1,numCycles), [], ...
                [], [], [],repmat(0,1,numCycles),[], ...
                @(x)numEntities(x,mu,type,d,N,METHOD,remoteP), options);
    
    centralFlow = zeros(1,numChains);
    cycleFlow = zeros(1,numCycles);     % this is now the same as DVs
    
    centralL = zeros(1,numChains);
    cycleL = zeros(1,numCycles);
    travelL = sum(x.*sum(d.*remoteP,2)'); %sum(x.*d);
    
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
        centralFlow(i) = sum(remoteP(:,i).*x');
        
        switch(METHOD)
          case 1,                       % M/M/1
                                        %disp('M/M/1 approximation');
            delta = centralFlow(i)/(mu(chains(i))-centralFlow(i));
            pk = 0;
            
          case 2,                       % M/M/1/N
            disp('FWR approximation');
            error('not implemented');
            rho = centralFlow(i)/mu(chains(i));
            pk = 0;
            if rho ~= 1
                % pk = ( rho^N*(1-rho) )/( 1 - rho^(N+1) );
                % Lq = rho/(1-rho) - rho*(N*rho^N + 1)/(1-rho^(N+1));
                w = (1 - (N+1)*rho^N + N*rho^(N+1))/(1-rho^(N+1));
                % see Whitt1984, eq. (2)
                delta = w*rho/(1-rho);
            else
                % pk = 1/(N+1);
                % Lq = N*(N+1)/( 2*(N+1) );
                delta = N/2;
            end
            %            delta = Lq + rho*(1-pk);
            
          case 3,                     % finite source
            disp('Finite source approximation');
            error('not implemented');
            rho = centralFlow(i)/mu(chains(i));
            a = zeros(1,N);
            for j = 1:N
                a(j) = nchoosek(N,j)*factorial(j)*rho^j;
            end
            p0 = 1/(1+sum(a));
            pk = p0*a(end);
            delta = p0*sum([1:N].*a);

          case 4,                       % D/M/1
            disp('D/M/1 approximation')
            error('not implemented');
            rho = centralFlow(i)/mu(chains(i));
            z = 0.5;
            oldZ = 0;
            while(abs(oldZ-z)>0.001)
                oldZ = z;
                z = exp(-(1/rho)*(1-z));
            end            
            delta = rho/(1-z);
            pk = 0;
        end
        
        centralL(i) = delta; %min(max(0,delta), N);
        centralW(i) = centralL(i)/( centralFlow(i)*(1-pk) );
    end
    
    for i = 1:numCycles
        % note that for finite-source, this is the *individual*
        % arrival rate, not the effective one...
        cycleFlow(i) = x(i); %sum(x(i:numCycles:end));

        switch(METHOD)
          case 1,                       % M/M/1
            delta = cycleFlow(i)/(mu(cycles(i))-cycleFlow(i));
            pk = 0;
            
          case 2,                       % M/M/1/N
            pk = 0;
            rho = cycleFlow(i)/mu(cycles(i));
            if rho ~= 1
                % pk = ( rho^N*(1-rho) )/( 1 - rho^(N+1) );
                % Lq = rho/(1-rho) - rho*(N*rho^N + 1)/(1-rho^(N+1));
                w = (1 - (N+1)*rho^N + N*rho^(N+1))/(1-rho^(N+1));
                % see Whitt1984, eq. (2)
                delta = w*rho/(1-rho);
            else
                % pk = 1/(N+1);
                % Lq = N*(N+1)/( 2*(N+1) );
                delta = N/2;
            end
            %            delta = Lq + rho*(1-pk);


          case 3,                       % finite-source
            rho = cycleFlow(i)/mu(cycles(i));
            a = zeros(1,N);
            for j = 1:N
                a(j) = nchoosek(N,j)*factorial(j)*rho^j;
            end
            p0 = 1/(1+sum(a));
            pk = p0*a(end);
            delta = p0*sum([1:N].*a);

          case 4,                       % D/M/1
            rho = cycleFlow(i)/mu(cycles(i));
            z = 0.5;
            oldZ = 0;
            while(abs(oldZ-z)>0.001)
                oldZ = z;
                z = exp(-(1/rho)*(1-z));
            end            
            delta = rho/(1-z);
            pk = 0;
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
        
      case 4,
        w = [centralW, cycleW];
    end
    q = [centralL, cycleL, travelL];
    
    %    disp(sprintf('System Throughput: %2.2f', sum(centralFlow)));

end


function T = throughput(lambda, weights)
    T = -sum(weights.*lambda);                   % negative to maximize!
end

function [C, Ceq] = numEntities(lambda, mu, type, d, N, METHOD, remoteP)

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
    travelL = sum(lambda.*sum(d.*remoteP,2)');

    for i=1:numChains
        centralFlow(i) = sum(remoteP(:,i).*lambda');

        switch(METHOD)
          case 1,                       % M/M/1
            delta = centralFlow(i)/(mu(chains(i))-centralFlow(i));
            pk = 0;
            
          case 2,                       % M/M/1/N
            error('not implemented');
            pk = 0;
            rho = centralFlow(i)/mu(chains(i));
            N;
            if rho ~= 1
                % pk = ( rho^N*(1-rho) )/( 1 - rho^(N+1) );
                % Lq = rho/(1-rho) - rho*(N*rho^N + 1)/(1-rho^(N+1));
                w = (1 - (N+1)*rho^N + N*rho^(N+1))/(1-rho^(N+1));
                % see Whitt1984, eq. (2)
                delta = w*rho/(1-rho);
            else
                % pk = 1/(N+1);
                % Lq = N*(N+1)/( 2*(N+1) );
                delta = N/2;
            end
            %            delta = Lq + rho*(1-pk);

          case 3,                       % finite-source
            error('not implemented');
            rho = centralFlow(i)/mu(chains(i));
            a = zeros(1,N);
            for j = 1:N
                a(j) = nchoosek(N,j)*factorial(j)*rho^j;
            end
            p0 = 1/(1+sum(a)); 
            delta = p0*sum([1:N].*a);

          case 4,                       % D/M/1
            error('not implemented');
            rho = centralFlow(i)/mu(chains(i));
            z = 0.5;
            oldZ = 0;
            while(abs(oldZ-z)>0.001)
                oldZ = z;
                z = exp(-(1/rho)*(1-z));
            end            
            delta = rho/(1-z);
            pk = 0;

        end
        centralL = centralL + delta; %min(max(0,delta), N);
    end
    

    for i = 1:numCycles
        cycleFlow(i) = lambda(i);%sum(lambda(i:numCycles:end));

        switch(METHOD)
          case 1,                       % M/M/1
            delta = cycleFlow(i)/(mu(cycles(i))-cycleFlow(i));
            pk = 0;
            
          case 2,                       % M/M/1/N
            rho = cycleFlow(i)/mu(cycles(i));
            pk = 0;
            if rho ~= 1
                % pk = ( rho^N*(1-rho) )/( 1 - rho^(N+1) );
                % Lq = rho/(1-rho) - rho*(N*rho^N + 1)/(1-rho^(N+1));
                w = (1 - (N+1)*rho^N + N*rho^(N+1) )/(1-rho^(N+1));
                % see Whitt1984, eq. (2)
                delta = w*rho/(1-rho);
            else
                % pk = 1/(N+1);
                % Lq = N*(N+1)/( 2*(N+1) );
                delta = N/2;
            end
            %            delta = Lq + rho*(1-pk);


          case 3,                       % finite-source
            rho = cycleFlow(i)/mu(cycles(i));
            a = zeros(1,N);
            for j = 1:N
                a(j) = nchoosek(N,j)*factorial(j)*rho^j;
            end
            p0 = 1/(1+sum(a));
            delta = p0*sum([1:N].*a);

          case 4,                       % D/M/1
            rho = cycleFlow(i)/mu(cycles(i));
            z = 0.5;
            oldZ = 0;
            while(abs(oldZ-z)>0.001)
                oldZ = z;
                z = exp(-(1/rho)*(1-z));
            end            
            delta = rho/(1-z);
            pk = 0;
        end
        
        cycleL = cycleL + delta; %min(max(0,delta), N);
    end
    
    if METHOD == 5                      % because DM1 can give
                                        % utilizations greater than 1
        C = [ centralL + cycleL + travelL - N];% total number only...
    
    else
        C = [ centralL + cycleL + travelL - N,... % total number
              centralFlow-mu(chains)+epsilon, ...  % stability of chains
              cycleFlow-mu(cycles)+epsilon];       % stability of cycles
    end
    
    Ceq = [];
end



