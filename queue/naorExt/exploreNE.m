function [NE, POA] = exploreNE(C,R, total, PLOT)
% NOTE: C = now a vector of size two

if nargin < 4
    PLOT = 0;
end

more off;

addpath '~/Documents/MATLAB/queueing/inst'

mu = [4 6];
Prob = [0 1; 1 0]; 

V = qnvisits(Prob);

%total = 1000;

i = 1;
U = zeros(2,total);
W = zeros(2,total);
Q = zeros(2,total);
X = zeros(2,total);

for N = 1:total
    [U(:,i), W(:,i), Q(:,i), X(:,i)] = qnclosed(N,1./mu,V);
    i = i + 1;
end

throughput = X(1,:);
wait = sum(W);
indivP1 = zeros(total, total);
indivP2 = zeros(total, total);

i = 1;
for n1 = 0:total
    j = 1;
    for n2 = 0:total
        idx = n1+n2;                    % total number in system;
                                        % index in performance metrics
        if idx <= total
            if n1 > 0 && n2 > 0
                indivP1(i,j) = n1/(n1+n2)*throughput(idx)*(R - C(1)* ...
                                                      wait(idx));
                indivP2(i,j) = n2/(n1+n2)*throughput(idx)*(R - C(2)* ...
                                                      wait(idx));
            elseif n1 == 0
                indivP1(i,j) = 0;
                if n2 == 0
                    indivP2(i,j) = 0;
                else
                    indivP2(i,j) = throughput(idx)*(R - C(2)* ...
                                                     wait(idx));
                end
            elseif n2 == 0
                indivP2(i,j) = 0;
                indivP1(i,j) = throughput(idx)*(R - C(2)* ...
                                                wait(idx));
            end

        else
            indivP1(i,j) = 0;
            indivP2(i,j) = 0;
        end
        j = j + 1;
    end
    i = i + 1;
end   

social = indivP1 + indivP2;

if PLOT
    mesh([0:total],[0:total],indivP1)
    axis([0 total/2 0 total/2]);
    title('Player 1');
    xlabel('n2');
    ylabel('n1');
    figure;
    contour([0:total],[0:total],indivP1)
    hold on;
    title('Player 1');
    xlabel('n2');
    ylabel('n1');
    [x,br1] = max(indivP1);
    br1 = br1 - 1;
    plot([0:total], br1, 'ks', 'MarkerSize', 10);
    axis([0 total/2 0 total/2]);

    figure;
    mesh([0:total],[0:total],indivP2)
    axis([0 total/2 0 total/2]);
    title('Player 2');
    xlabel('n2');
    ylabel('n1');
    figure;
    contour([0:total],[0:total],indivP2)
    hold on;
    title('Player 2');
    xlabel('n2');
    ylabel('n1');
    [x,br2] = max(indivP2');
    br2 = br2 - 1;
    plot(br2, [0:total], 'ko', 'MarkerSize', 15);

    axis([0 total/2 0 total/2]);

    figure;
    mesh([0:total],[0:total],social)
    axis([0 total/2 0 total/2]);
    title('Total');
    xlabel('n2');
    ylabel('n1');
    
    figure;
    contour([0:total],[0:total],social)
    hold on;
    plot([0:total], br1, 'ks', 'MarkerSize', 10);
    plot(br2, [0:total], 'ko', 'MarkerSize', 15);
    title('Total');
    xlabel('n2');
    ylabel('n1');
    axis([0 total/2 0 total/2]);

    [tmp, max1] = max(social);
    [tmp, max2] = max(max(social));
    plot(max2-1, max1(max2)-1, 'ro', 'MarkerSize', 20);

end

[x,br1] = max(indivP1);
br1 = br1 - 1;
[x,br2] = max(indivP2');
br2 = br2 - 1;
[tmp, max1] = max(social);
[tmp, max2] = max(max(social));

NEcount = 0;
for i =0:total
    if i == br1(br2(i+1)+1)
        NEcount = NEcount + 1;
        NE1 = br1(br2(i+1)+1);
        NE2 = br2(i+1);
    end
end

if NEcount == 0
    NE1 = 0;
    NE2 = 0;
elseif PLOT
    plot(NE2, NE1, 'r*', 'MarkerSize', 15);
end

% OUTPUTS
POA = social(max2,max1(max2))/social(NE1+1,NE2+1);
NE = NEcount;

disp(sprintf('Optimal: %f \t Nash(%d): %f \t PoA: %f',social(max2, ...
                                                  max1(max2)), ...
             NEcount, social(NE1+1,NE2+1), POA));