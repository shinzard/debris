% This function performs disrete-event queueing network simulation
% for single-chain closed queueing networks using the same
% interface as 'qnclosed' in the queueing toolbox:
%   FUNCTION [U,W,Q,X] = QUEUESIM(N,S,P,T,WARMUP)
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
%     T - [optional - def. 5000] simulation end time
%     WARMUP - [optional - def. 1000] warmup period
%
% 2 Nov 2012
% J.Brooks
function [u,w,q,x] = queueSim(N, S, P, T, warmup)

    forgetFactor = 0.9;
    
    if nargin < 5
        T = 5000;                      % same units as S
        warmup = 1000;                 % same units as S
    end
    
    if sum(P) ~= 1 && P ~= 0 && P ~= -1 && P ~= -2
        error('Probability vector error');
    end

    % simulation variables
    eventList = NaN*zeros(1,sum(N));    % next time for each entity
    t = 0;                              % current simulation time
    lastT = 0;                          % last sim. time

    % initialization of server metrics
    busy = zeros(1,length(S));          % busy counter for each
                                        % server
    numServed = zeros(1,length(S));
    cumWaitTime = zeros(1,length(S));
    cumLength = zeros(1,length(S));
    currLengths = zeros(1,length(S));
    waitFilt = S;                       % initialize with expected
                                        % service times
    
    % initialization of entity variables
    class = zeros(1,sum(N));
    enterTime = zeros(1,sum(N));
    
    if P == 0                           % deterministic routing
        tmp = 1;
        classType = 1;

        for i = 1:length(N)
            class(tmp:tmp + N(i)-1) = repmat(classType, 1, N(i));
            classType = classType + 1;
            tmp = tmp + N(i);
        end
    end

    P = cumsum(P);                      % to use with unif random number

    location = ones(1,sum(N));          % start all entities at
                                        % first server

    position = [1:sum(N)];              % entitity in position one
                                        % is at server

    eventList(1) = exprnd(S(1));        % end of first service
    
    
    while t < T
        for i = 1:length(S)
            currLengths (i)= length(find(location == i));
        end

        lastT = t;
        [t,idx] = min(eventList);       % get next event (end of
                                        % service time for entity
                                        % idx)
        eventList(idx) = NaN;

        % address leaving queue
        queue = location(idx);
        neighbors = find(location == queue);

        if t > warmup
            cumLength = cumLength + currLengths*(t-lastT);
        end

        if ~isempty(neighbors)
            position(neighbors) = max(0, position(neighbors) - 1); 
            nextUp = neighbors(find(position(neighbors)==1));
        
            if ~isempty(nextUp)         % should always be non-empty...
                nextService = exprnd(S(queue));
                eventList(nextUp) = t + nextService; 
            
                % Add up server metrics
                if t > warmup
                    busy(queue) = busy(queue) + nextService;
                    numServed(queue) = numServed(queue) + 1;
                    tmpWait = (t - enterTime(nextUp) + nextService);
                    cumWaitTime(queue) = cumWaitTime(queue) + tmpWait;
                    waitFilt(queue) = (1-forgetFactor)*waitFilt(queue) ...
                        + forgetFactor*tmpWait;
                end
            end
        end
        
        % Address Entering queue (and routing)
        if queue == 1
            if class(idx) == 0          % probabilistic or rational
                                        % routing
                if P == -1              % rational
                                        %[tmp, nextQueue] =
                                        %min(waitFilt(2:end));
                    [tmp, nextQueue] = min(currLengths(2:end).*S(2: ...
                                                                 end));
                    nextQueue = nextQueue + 1;

                elseif P == -2          % heuristic routing
                    [tmp, nextQueue] = min(currLengths(2:end).*(S(2: ...
                                                                 end).^(0.5)));
                    nextQueue = nextQueue + 1;
                    
                    
                else                    % probabilistic
                    tmp = unifrnd(0,1);
                    nextQueue = find(tmp < P, 1, 'first')+1;
                end

            else                        % deterministic routing
                nextQueue = class(idx) + 1;% these are offset because
                                           % the central server is '1'
            end

        else
            nextQueue = 1;              % always loop back
        end
        
        enterTime(idx) = t;

        neighbors = find(location == nextQueue);
        location(idx) = nextQueue;

        if isempty(neighbors)           % Start service right away
            position(idx) = 1;
            nextService = exprnd(S(nextQueue));
            eventList(idx) = t + nextService;
            
            % Add up server metrics
            if t > warmup
                busy(nextQueue) = busy(nextQueue) + nextService;
                numServed(nextQueue) = numServed(nextQueue) + 1;
                tmpWait = (t - enterTime(idx) + nextService);
                cumWaitTime(nextQueue) = cumWaitTime(nextQueue) + tmpWait;
                waitFilt(nextQueue) = (1-forgetFactor)*waitFilt(nextQueue) ...
                    + forgetFactor*tmpWait;
            end
        else                            % get in line
            position(idx) = max(position(neighbors)) + 1;
        end
    end
    
    % Calculate performance metrics for run;
    u = busy./(T-warmup);
    w = cumWaitTime./numServed;
    q = cumLength./(T-warmup);
    x = u./S;
    
    numServed
