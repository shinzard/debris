close all;
clear all;
more off;

mu1 = 10;
mu2V = [1:0.2:25];
mu3V = [5:5:40];
NV = 24; %[6:3:24];

% System configuration
P = zeros(2,3,2,3);
P(1,1,1,3) = 1;
P(1,3,1,1) = 1;
P(2,2,2,3) = 1;
P(2,3,2,2) = 1;
V = qnvisits(P);

% Results vectors
NUM_TESTS = length(mu2V)*length(mu3V)*length(NV);
u = zeros(NUM_TESTS,3);
w = zeros(NUM_TESTS,3);
q = zeros(NUM_TESTS,3);
x = zeros(NUM_TESTS,3);
n1 = zeros(NUM_TESTS,1);
mu3r = zeros(NUM_TESTS,1);
mu2r = zeros(NUM_TESTS,1);
nr = zeros(NUM_TESTS,1);

test = 1;

for mu3 = mu3V
    disp(sprintf('Running mu3 = %d', mu3));
    for mu2 = mu2V
        disp(sprintf('Running mu2 = %d', mu2));
        for N = NV
            opt = 0;
            mu3r(test) = mu3; mu2r(test) = mu2; nr(test) = N;
            S = repmat([1/mu1, 1/mu2, 1/mu3],2,1);

            % choose best among these
            for N1 = [1:N-1]
                [ut,wt,qt,xt] = qnclosed([N1,N-N1],S,V);
                if sum(xt(:,3))>opt
                    u(test,:) = sum(ut); w(test,:) = sum(wt); q(test,:) = sum(qt);
                    x(test,:) = sum(xt); n1(test) = N1;
                    opt = sum(xt(:,3));
                end
            end
            disp(sprintf('Optimal: [%d,%d]', n1(test), N-n1(test)));
            test = test + 1;
        end
    end
end

for N = NV
    idx = find(nr == N);
    figure, title(sprintf('Optimal Ratio; N=%d',N));
    hold on;
    figRatio = gcf;
    figure, title(sprintf('Throughput; N=%d',N));
    hold on;
    figThr = gcf;
    for mu3 = mu3V
        idx3 = find(mu3r(idx) == mu3);
        figure(figRatio), plot(10./mu2r(idx(idx3)), n1(idx(idx3))./ ...
                               nr(idx(idx3)), 'b');
        %        text(5,0.8, sprintf('mu3=%d', mu3));
        figure(figThr), plot(10./mu2r(idx(idx3)), x(idx(idx3),3), ...
                             'b');
        %        endtext(5,3, sprintf('mu3=%d', mu3));
    
    end
end        