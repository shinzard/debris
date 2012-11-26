close all;
clear all;
more off;

mu1 = 10;
mu2V = [5,10,15]; %[1:2:25];
mu3V = 15; %[5:5:40];
NV = 12; %[6:3:24];

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
plotIdx = 1;

figure(1);
subplot(4,1,1), title('Utilizations'), hold on;
subplot(4,1,2), title('Wait Times'), hold on;
subplot(4,1,3), title('Number of Customers'), hold on;
subplot(4,1,4), title('Throughput'), hold on;

for mu3 = mu3V
    disp(sprintf('Running mu3 = %d', mu3));
    for mu2 = mu2V
        disp(sprintf('Running mu2 = %d', mu2));
        for N = NV
            opt = 0;
            mu3r(test) = mu3; mu2r(test) = mu2; nr(test) = N;
            S = [1/mu3, 1/mu1, 1/mu2];
            S2 = repmat([1/mu1, 1/mu2, 1/mu3],2,1);

            % choose best among these
            for N1 = [1:N-1]
                [ut2,wt2,qt2,xt2] = qnclosed([N1,N-N1],S2,V);
                [ut,wt,qt,xt] = queueSim([N1,N-N1],S,0,3000,500);

                figure(1), subplot(4,1,1);
                plot(repmat(plotIdx,1,3), sum(ut2), 'b*', 'markersize', ...
                     15);
                hold on;
                plot(repmat(plotIdx,1,3), ut, 'r^', 'markersize', 10);

                subplot(4,1,2);
                plot(repmat(plotIdx,1,3), [wt2(1,1), wt2(2,2), (N1*ut2(1,3) + ...
                                                          (N-N1)* ...
                                                          wt2(2,3))/N], ...
                     'b*', 'markersize', 15);
                
                plot(repmat(plotIdx,1,3), wt, 'r^', 'markersize', 10);

                subplot(4,1,3), plot(repmat(plotIdx,1,3), sum(qt2), 'b*', 'markersize', 15);
                plot(repmat(plotIdx,1,3), qt, 'r^', 'markersize', 10);

                subplot(4,1,4), plot(repmat(plotIdx,1,3), sum(xt2), 'b*', 'markersize', 15);
                plot(repmat(plotIdx,1,3), xt, 'r^', 'markersize', 10);
                plotIdx = plotIdx + 1;
                
                if sum(xt(1))>opt
                    u(test,:) = ut; w(test,:) = wt; q(test,:) = qt;
                    x(test,:) = xt; n1(test) = N1;
                    opt = xt(1);        % server 1 is central by
                                        % definition (note
                                        % difference from earlier
                                        % MVA analysis)
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