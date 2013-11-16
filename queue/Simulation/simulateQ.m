% This function performs disrete-event queueing network simulation
% for general bipartite closed queueing networks using the same
% interface as 'qnclosed' in the queueing toolbox:
%   FUNCTION [U,W,Q,X] = QUEUESIM(N,S,P,T,WARMUP)
%   Outputs:
%     U - utilization vector
%     W - response time (W = Wq + S)
%     Q - average number of customers (L = Lq + p)
%     X - throughput
%
%   Inputs:
%     N - initial vector of entities by routing class (length = # parallel cycles)
%     S - vector of mean service times
%     P - routing probability vector (0 if determinstic routing by
%     class, -1 for dispatcher routing)
%     (DISPATCH/COMP)PARAMS - parameters to be passed on
%     REMOTEP - probabilistic routing from remote to central sites
%     (only for P == -1, 0)
%     T - [optional - def. 5000] simulation end time
%     WARMUP - [optional - def. 1000] warmup period
%
% 4 Nov 2013
% J.Brooks
% Last Modified: 12 Nov. 2013

function [u,w,q,x] = simulateQ(N, S, type, P, dij, dispatchParams, compParams, remoteP, ...
                               T, warmup)

    if nargin < 7
        T = 5000;                      % same units as S
        warmup = 1000;                 % same units as S
    end
    
    DEBUG = 0;
    
    % simulation variables
    eventList = NaN*zeros(1,sum(N) + 1);% next time for each entity
                                        % + dispatcher
    t = 0;                              % current simulation time
    lastT = 0;                          % last sim. time
    cumHist = zeros(sum(N),sum(N));
    today = zeros(sum(N),sum(N));

    % initialization of server metrics
    numStations = size(S,2);
    numCentral = length(find(type==1));
    busy = zeros(1,numStations);        % busy counter for each
                                        % server

    % Team attribute vectors
    teamFam = zeros(1,numStations-numCentral);
    teamSize = zeros(1,numStations-numCentral);


    numServed = zeros(1,numStations);
    cumWaitTime = zeros(1,numStations);
    cumLength = zeros(1,numStations);
    currLengths = zeros(1,numStations);
    waitFilt = S;                       % initialize with expected
                                        % service times
    
    % initialization of entity variables
    class = ones(1,sum(N));
    enterTime = zeros(1,sum(N));
    
    if length(N)>1 %P == 0                           % deterministic routing
        tmp = 1;
        classType = 1;

        for i = 1:length(N)
            class(tmp:tmp + N(i)-1) = repmat(classType, 1, N(i));
            classType = classType + 1;
            tmp = tmp + N(i);
        end
        
        if P <= 0                       % create deterministic
                                        % routing matrix
            P = zeros(length(N),numStations,length(N),numStations);
            for i = 1:length(N)
                for j = 1:numCentral
                    P(i,j,i,i+numCentral) = 1;   % always return to remote site
                    P(i,i+numCentral,i,j) = remoteP(i,j); % choose central
                                                 % site based on
                                                 % input probability
                end
            end
            
            %            reshape(P(1,:,1,:),8,8)


        end
    end

    %P = cumsum(P);                      % to use with unif random number
    location = zeros(1,sum(N));

    if length(N)>1
        for i = 1:length(N)             % for each class, find
                                        % feasible starting point
            idx = find(class == i);
            tmp = mod(find(P(i,:,i,:),1,'first'),size(S,2));
            if tmp > 0
                location(idx) = tmp;
            else
                location(idx) = size(S,2);
            end
        end
    else
        location = location + 1;        % start all at first node
    end

    % Assign initial positions (1 = at server)
    positionTmp = ones(1,size(S,2));
    for i = 1:length(location)
        position(i) = positionTmp(location(i));
        positionTmp(location(i)) = positionTmp(location(i)) + 1;
    end

    % Start service
    idx = find(position == 1);
    for i = 1:length(idx)
        eventList(idx(i)) = exprnd(1/S(location(idx(i))));
    end
    
    % First dispatcher decision
    eventList(end) = exprnd(dispatchParams.inter);
    
    % ----------------------------------------
    % RUN SIMULATION
    % ----------------------------------------
    while t < T
        %        class
        %        eventList
        %        position
        %        currLengths
        % Update familiarity
        [teamFam,teamSize,cumHist,today] = familiarity(t, lastT, cumHist, today, ...
                                            class,teamFam);
        
        for i = 1:numStations
            currLengths (i)= length(find(location == i));
        end

        lastT = t;
        [t,idx] = min(eventList);       % get next event (end of
                                        % service time for entity
                                        % idx)
        eventList(idx) = NaN;
        
        if idx == sum(N)+1              % dispatcher event
            class = dispatcher(class, dispatchParams, waitFilt); % reallocate

            % Set next dispatcher decision
            eventList(end) = t + exprnd(1/dispatchParams.inter);
            continue;
        end
        
        % For debuggins....
        if idx == 2 && DEBUG
            disp(sprintf('t = %2.2f; location: %d; pos: %d, class: %d',t, location(idx), ...
                         position(idx), class(idx)));
        end

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
                if queue > numCentral
                    nextService = ...
                        exprnd(1/compositionEffect(teamFam(queue- ...
                               numCentral),teamSize(queue-numCentral),S(queue),compParams));
                else
                    nextService = exprnd(1/S(queue));
                end
                eventList(nextUp) = t + nextService; 
            
                % Add up server metrics
                if t > warmup
                    busy(queue) = busy(queue) + nextService;
                    numServed(queue) = numServed(queue) + 1;
                    tmpWait = (t - enterTime(nextUp) + nextService);
                    cumWaitTime(queue) = cumWaitTime(queue) + tmpWait;
                    waitFilt(queue) = (1-dispatchParams.forgetFactor)*waitFilt(queue) ...
                                        + dispatchParams.forgetFactor*tmpWait;
                end
            end
        end
        
        % Address Entering queue (and routing)
        tmp = unifrnd(0,1);
        if length(size(P))>2    % many entity classes
                                % (pure, mixed)
            nextP = P(class(idx),queue,class(idx),:);
        else                    % only one class (strict prob)
            nextP = P(queue,:);
        end

        nextQueue = find(tmp < cumsum(nextP), 1, 'first');
        % disp(sprintf('Class: %d; Leaving %d, entering %d\n', ...
        %              class(idx),queue,nextQueue));
        % input('Next');
            
        enterTime(idx) = t;
        
        neighbors = find(location == nextQueue);
        location(idx) = nextQueue;

        if isempty(neighbors)           % Start service right away
            position(idx) = 1;
            if nextQueue > numCentral   % remote site (with team effects)
                nextService = ...
                    exprnd(1/compositionEffect(teamFam(nextQueue-numCentral),teamSize(nextQueue-numCentral),S(nextQueue),compParams));                
            else                        % central site
                nextService = exprnd(1/S(nextQueue));
            end

            eventList(idx) = t + nextService;
            
            % Add up server metrics
            if t > warmup
                busy(nextQueue) = busy(nextQueue) + nextService;
                numServed(nextQueue) = numServed(nextQueue) + 1;
                tmpWait = (t - enterTime(idx) + nextService);
                cumWaitTime(nextQueue) = cumWaitTime(nextQueue) + tmpWait;
                % waitFilt(nextQueue) = (1-forgetFactor)*waitFilt(nextQueue) ...
                %                       + forgetFactor*tmpWait;
            end
        else                            % get in line
            position(idx) = max(position(neighbors)) + 1;
        end
        
        %        figure(1);
        %        plot(t, teamSize(1), 'b.')
        %        hold on;
    end
    
    % Calculate performance metrics for run;
    u = min(1,busy./(T-warmup));
    w = cumWaitTime./numServed;
    q = cumLength./(T-warmup);
    x = numServed./(T-warmup);%u./S(1,:);