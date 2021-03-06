function [Nvec, p] = routing(type, mu, d, f , N, METHOD)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    chains = find(type == 1);
    cycles = find(type == 2);
    numChains = length(chains);
    numCycles= length(cycles);
    decisionVars = length(f);
    numStations = length(mu);
    centralFlow = zeros(1,numChains);
    cycleFlow = zeros(1,numCycles);
    centralL = zeros(1,numChains); 
    cycleL = zeros(1,numCycles); 

    travelL = f.*d;

    for i=1:numChains
        % note that for finite-source, this is the *individual*
        % arrival rate, not the effective one...
        centralFlow(i) = sum(f((i-1)*numCycles + 1:i* ...
                               numCycles));
        switch(METHOD)
          case 1,                       % M/M/1
            disp('M/M/1 approximation');
            delta = centralFlow(i)/(mu(chains(i))-centralFlow(i));
            pk = 0;
            
          case 2,                       % M/M/1/N
            disp('FWR approximation');
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
        cycleFlow(i) = sum(f(i:numCycles:end));
        
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

    % for j=1:numChains
    %     centralFlow(j) = sum(f((j-1)*numCycles + 1:j* ...
    %                            numCycles));
    %     switch(METHOD)
    %       case 1,                       % M/M/1
    %         delta = centralFlow(j)/(mus(chains(j))-centralFlow(j));
    %         pk = 0;
            
    %       case 2,                       % M/M/1/N
    %         rho = centralFlow(j)/mus(chains(j));
    %         if rho < 1
    %             pk = ( rho^N*(1-rho) )/( 1 - rho^(N+1) );
    %             Lq = rho/(1-rho) - rho*(N*rho^N + 1)/(1-rho^(N+1));
    %         else
    %             pk = 1/(N+1);
    %             Lq = N*(N+1)/( 2*(N+1) );
    %         end
    %         delta = Lq + rho*(1-pk);
          
    %       case 3,                       % finite-source
    %         rho = centralFlow(j)/mus(chains(j));
    %         a = zeros(1,N);
    %         for k = 1:N
    %             a(k) = nchoosek(N,k)*factorial(k)*rho^k;
    %         end
    %         p0 = 1/(1+sum(a));
    %         delta = p0*sum([1:N].*a);
    %     end
    %     centralL(j) = delta; %min(max(0,delta), N);
    % end
    
    % for j = 1:numCycles
    %     cycleFlow(j) = sum(f(j:numCycles:end));

    %     switch(METHOD)
    %       case 1,                       % M/M/1
    %         delta = cycleFlow(j)/(mus(cycles(j))-cycleFlow(j));
    %         pk = 0;
            
    %       case 2,                       % M/M/1/N
    %         rho = cycleFlow(j)/mus(cycles(j));
    %         if rho < 1
    %             pk = ( rho^N*(1-rho) )/( 1 - rho^(N+1) );
    %             Lq = rho/(1-rho) - rho*(N*rho^N + 1)/(1-rho^(N+1));
    %         else
    %             pk = 1/(N+1);
    %             Lq = N*(N+1)/( 2*(N+1) );
    %         end
    %         delta = Lq + rho*(1-pk);
    %       case 3,                       % finite-source
    %         rho = cycleFlow(j)/mus(cycles(j));
    %         a = zeros(1,N);
    %         for k = 1:N
    %             a(k) = nchoosek(N,k)*factorial(k)*rho^k;
    %         end
    %         p0 = 1/(1+sum(a));
    %         delta = p0*sum([1:N].*a);
    %     end
        
    %     cycleL(j) = delta; %min(max(0,delta), N);
    % end

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

    while remaining > length(Nvec)      %  only happens with FWR
                                        %  method when at least one
                                        %  server has a utilization
                                        %  of 1
        NvecOrig
        Nvec
        remaining

        warning('invoking strange routing modification');
        Nvec = Nvec + 1;
        [x,idx] = sort(NvecOrig - Nvec, 'descend');
        remaining = N - sum(Nvec);

    end

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

