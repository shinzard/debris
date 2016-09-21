% This function defines the dispatcher decision logic. It is called
% from simulateQ at random intervals as specified based on current
% system state and the parameters.
% 
% 4 Nov 2013
% J.Brooks
%  
%  Last Modified: 2 March 2014
%  
function newAssign = dispatcher(oldAssign, location, param, prices, ...
                                state, numCentral, remoteP, mus, ...
                                type, dij, N, history)
n = size(state,1);
newAssign = oldAssign;

switch(param.method)
  case 'Fixed'                          % FIXED
    newAssign = oldAssign;
  case 'Bounded'                        % Bounded Rationality
    ak = zeros(1,n*(n-1));              % all possible switches
    idx = 1;
    from = zeros(1,n*(n-1));
    to = zeros(1,n*(n-1));

    for i = 1:size(state,1)
        for j = 1:size(state,1)
            if i ~= j                   % avoid the no-change case
                ak(idx) = log(( prices(i)/sum(state(i,:).*remoteP(i,:)) )/ ...
                              ( prices(j)/sum(state(j,:).*remoteP(j,:)) ));
                to(idx) = i;
                from(idx) = j;
                idx = idx + 1;
            end
        end
    end
    
    p = exp(ak/param.T)./sum(exp(ak/param.T));
    % DEBUG
    %    figure, stem(p);
    tmp = unifrnd(0,1);
    dec = find(tmp < cumsum(p), 1, 'first');
    trucksToMove = find(oldAssign == from(dec));

    % DEBUGGING
    %     figure, plot(p, 'b.'), hold on;
    %     plot(dec, p(dec), 'rs');
    %     disp(sprintf('Dispatch: %d-->%d\t ak: %2.2f', from(dec), to(dec), ak(dec)));
    
    if isempty(trucksToMove)            % nothing to do...
        return;
    end
    
    % look for trucks at central site (to match ARENA model...)
    atCentral = trucksToMove(find(location(trucksToMove) <= numCentral));
    
    if ~isempty(atCentral)
        newAssign(randsample(atCentral,1)) = to(dec);
    else
        newAssign(randsample(trucksToMove, 1)) = to(dec);
    end
    
  case 'Optimal'

    [ut,wt,qt,xt,et]=optimalAssignmentFixedProb(mus, type, dij, ...
                                                repmat(1/6,1,6), ...
                                                sum(N), 1, remoteP); 
    
    Nvec = routingFixed(type, mus, dij, xt, sum(N), 1, remoteP);

    x0 = zeros(6*24, 1);

    for i = 1:24
        x0(oldAssign(i)+6*(i-1)) = 1;
    end
    
    x = assignment(history, Nvec, x0);

    newAssign = find(abs(reshape(x,6, 24)-1)<0.01)-[0:6:143]';
    
    % newAssign =  [repmat(1,1,4), repmat(2,1,5), repmat(3,1,5), ...
    %               repmat(4,1,4), repmat(5,1,5), repmat(6,1,1)];
end
