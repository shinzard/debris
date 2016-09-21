% This function determines the maximum total working history
% assignment given required loop counts
%
% 18 Apr 2013
% J.Brooks

function x = assignment(history, optFlows, x0)

numTrucks = size(history,1);
numTeams = length(optFlows);

if numTrucks ~= sum(optFlows)
    error('size mismatch');
end

numVars = numTrucks*numTeams;
x = zeros(1,numVars);
H = zeros(numVars,numVars);

%------------------------------
% Build H matrix (min --> neg)
%------------------------------
for i = 1:numTrucks-1
    hTmp = [];
    histDiag = diag(history,i);
    for j = 1:length(histDiag)
        hTmp = [hTmp, -repmat(histDiag(j),1,numTeams)];
    end
    H = H + diag(hTmp, i*numTeams);
end

%------------------------------
% Build equality constraints
%------------------------------
Aeq = zeros(numTrucks+numTeams, numVars);

% Assigment constraints
for i = 1:numTrucks
    Aeq(i,:) = [repmat(zeros(1,numTeams), 1, i-1), ones(1,numTeams), ...
                repmat(zeros(1,numTeams), 1, numTrucks - i)];
end

% Team constraints
for i = 1:numTeams
   Aeq(i + numTrucks,:) = repmat([repmat(0, 1, i-1), ...
                     1, repmat(0, 1, numTeams - i)], ...
                                1, numTrucks);
end

%------------------------------
% Setup problem and solve
%------------------------------
F = zeros(1,numVars);
Beq = [ones(1,numTrucks),optFlows]';
LB = zeros(1,numVars);
UB = ones(1,numVars);

% Need to replace with call to CPLEX...
%x = quadprog(H,F,[],[],Aeq,Beq,LB,UB);

% CPLEX notes: need to ensure H is symmetric and psd...
delta = -min(min(H));
H = 1/2*(H+H');% + delta*eye(numVars)
opts = cplexoptimset('Display', 'iter', 'Diagnostics', 'off', 'MaxIter', ...
                     3000, 'MaxNodes', 5000);

[x,fval,flag,output] = cplexmiqp(H,F,[],[],Aeq,Beq,[],[],[],[],[], ...
                                 repmat('B',1,numVars),x0, opts);
%H

%[x,fval,flag,output] = cplexmiqp(H,F,[],[],Aeq,Beq,[],[],[],LB,UB, ...
%                                 [],x0, opts);