addpath('../');
clear all;
close all;

minFlowI = 1e6;
minFlowIw = [];

% Kobayashi and Gerla, 1983 (Example 1)
mus = [4 2 1 0.5];
type = [1 2 2 2];

NUM_ENTITIES = 5;

N = 101;
w1 = linspace(0,1,N);
w2 = linspace(0,1,N);
[w1, w2] = ndgrid(w1,w2);
w3 = 1 - w1 - w2;
bad = (w1+w2 > 1); w1(bad) = NaN; w2(bad) = NaN; w3(bad) = NaN;

for i = 1:N
    for j = 1:N
        if ~isnan(w1(i,j))
            [u,w,q,x,e]=optimalAssignmentWeights(mus,type,0,[w1(i,j),w2(i,j),w3(i,j)],NUM_ENTITIES,1);
            %[Nvec, p] = routing(type, mus, 0, x, NUM_ENTITIES, 1);
        
            throughput(i,j) = sum(x);
            
            % Flows
            inequity(i,j) = sqrt(sum(((x./mean(x)-1).^2)));
            
            if inequity(i,j) < minFlowI
                minFlowIw = [w1(i,j),w2(i,j),w3(i,j)];
                minFlowI = inequity(i,j);
            end
            
            % Wait
            inequity2(i,j) = sqrt(sum(((w(2:end)./mean(w(2:end))-1).^2)));
            
            % Length/Utilization
            inequity3(i,j) = sqrt(sum(((q(2:end-1)./mean(q(2:end-1))-1).^2)));
            
            % Utilization
            inequity4(i,j) = sqrt(sum(((u(2:end)./mean(u(2:end))-1).^2)));
        else
            throughput(i,j) = NaN;
            inequity(i,j) = NaN;
            inequity2(i,j) = NaN;
            inequity3(i,j) = NaN;
            inequity4(i,j) = NaN;
        end
    end
end

%save(['simplexResults-',num2str(NUM_ENTITIES),'_', datestr(now)]);

%simplexPlot;
