close all;
clear all;
more off;

type = [1, 1, 2, 2, 2];
central = find(type == 1);
parallel = find(type == 2);

filenameroot = '../2-3networks'

configs = dlmread([filenameroot, 'config.txt']);

NUM_TESTS = configs(1,1);
configs = configs(2:end,:);             % thow out first row

% Results vectors
%waitEq = zeros(1,NUM_TESTS);
%lengthEq = zeros(1,NUM_TESTS);
flowEq = zeros(1,NUM_TESTS);
opt = zeros(1,NUM_TESTS);
optTOconvex = zeros(1,NUM_TESTS);
optTOiter = zeros(1,NUM_TESTS);
%Imu = zeros(1,NUM_TESTS);
Iopt = zeros(1,NUM_TESTS);
ItoConvex = zeros(1,NUM_TESTS);
ItoIter = zeros(1,NUM_TESTS);
Pmatrix = zeros(NUM_TESTS,6);

for i = 1:NUM_TESTS
    disp(sprintf('Config: %d', i));
    
    mus = configs(i,1:5);
    d = configs(i,6:end-1);
    d = reshape(d, length(parallel), length(central));
    N = configs(i,end);

    %    Imu(i) = sqrt(sum((mus(parallel)./mean(mus(parallel))-1).^2))/sqrt(3*2);
    
    % Generate random fixed routing probabilities...
    pr = rand(length(parallel),1);
    pr = [pr, 1-pr];
    Pmatrix(i,:) = reshape(pr,1,6);

    % Optimal
    weightsEff = repmat(1/length(parallel),1,length(parallel)); % maximum
                                                                % effectiveness
    
    [u,w,q,x,exit]=optimalAssignmentFixedProb(mus,type,d,weightsEff,N, ...
                                              1,pr);

    opt(i) = sum(x);

    Iopt(i) = sqrt(sum((x./mean(x)-1).^2)/(3*2));
    
    % Wait-Equity
    %    waitEq(i) = equityLinesearch2(type, mus, d, pr, 0, N);
    
    % Length-Equity
    %    lengthEq(i) = equityLinesearch2(type, mus, d, pr, ...
    %                                       mus(parallel)/sum(mus(parallel)), N);
    
    
    % Flow-Equity
    flowEq(i) = equityLinesearch2(type, mus, d, pr, repmat(1/length(parallel), ...
                                                      1, length(parallel)), N);

    parallelFlows = flowEq(i).*repmat(1/length(parallel), 1, ...
                                      length(parallel));
    
    centralFlows = parallelFlows*pr;
    
    % --------------------------------------------------
    % COMPARE TRADEOFF METHODS
    % --------------------------------------------------

    weightsEq = [sum(pr(1,:).*mus(central)./ ...
                     (mus(central)-centralFlows).^2), ...
                 sum(pr(2,:).*mus(central)./ ...
                     (mus(central)-centralFlows).^2), ...
                 sum(pr(3,:).*mus(central)./ ...
                     (mus(central)-centralFlows).^2)] + ...
        mus(parallel)./(mus(parallel)-parallelFlows).^2 + ...
        sum(pr.*d,2)';

    [u,w,q,x,exit]=optimalAssignmentFixedProb(mus,type,d,weightsEff*0.5 ...
                                              + weightsEq*0.5,N, ...
                                              1,pr);

    optTOconvex(i) = sum(x);

    ItoConvex(i) = sqrt(sum((x./mean(x)-1).^2)/(3*2));
    
    w1 = weightsEff;

    minI = 1e6;
    
    for j = 1:100
        [u,w,q,x,e]=optimalAssignmentFixedProb(mus, type, d, w1, N, ...
                                               1, pr);
        Itmp = sqrt((x./mean(x) - 1).^2);
        n = (x<mean(x))*2-1;
        w1 = w1 + 0.01.*Itmp.*n;
        w1 = w1/sum(w1);

        % take minimum solution from the 20...
        if sqrt(sum((x./mean(x)-1).^2)/(3*2)) < minI
                optTOiter(i) = sum(x);
                ItoIter(i) = sqrt(sum((x./mean(x)-1).^2)/(3*2));
                minI = sqrt(sum((x./mean(x)-1).^2)/(3*2));
        end
    end



end

save(['equityLineSearch2x3-',num2str(N),'_', datestr(now)]);

plotEqCompare;