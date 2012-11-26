close all;
clear all;

%lambda = 1;
mu1 = 18;
mu2 = 25;                               % picked to allow up to
                                        % p=0.8 without instability
%%%
%figure(1);
%for p = 0:0.01:0.95
%    L1 = lambda/((1-p)*mu1 - lambda);
%    L2 = p*lambda/((1-p)*mu2 - p*lambda);
%    L = L1 + L2;
%    plot(p, L, 'b*', 'markersize', 20);
%    hold on;
%    plot(p, L1, 'r*', 'markersize', 15);
%    plot(p, L2, 'g*', 'markersize', 15);
%end

P = [0,1; 1,0];
V = qnvisits(P);
S = [1/mu1, 1/mu2];

% For new approx method
l = 0.001;

for N = 1:30;
    % New Approx Method
    tic;                                % can we improve with lower
                                        % lambda (higher p)??
    
    %% Assumes lambda = 1:
    % p = min(roots([-mu1-N*(mu1*mu2 + mu1),mu1-mu2-2-N*(1+mu2-mu1-2*mu1*mu2), ...
    %                   mu2-N*(mu1*mu2-mu2)]));
    
    pRoot = roots([-l*(mu1+l)-N*((mu1+l)*mu2 + l*(mu1+l)),l*((mu1+l)-mu2-2*l)-N*(l^2+l*mu2-l*(mu1+l)-2*(mu1+l)*mu2), ...
                   l*mu2-N*((mu1+l)*mu2-l*mu2)]);
    p = min(pRoot);
    w1 = [(1-p)/((1-p)*(mu1+l)-l), (1-p)/((1-p)*mu2-p*l)];
    x1 = p*l/(1-p);
    t1 = toc;

    % check
    %    l1 = l/((1-p)*mu1 - l); 
    %    l2 = p*l/((1-p)*mu2 - p*l);
    %    l1+l2

    % MVA
    tic;
    [u,w,q,x] = qnclosed(N,S,V);
    t2 = toc;
    
    % Old Approx Method
    tic;
    x2 = min(roots([1, 1-N-mu1-mu2, mu1*mu2+mu2+N*(mu1+mu2), ...
                -N*mu1*mu2]));
    if (!isreal(x2))
        disp(sprintf('Non-real old approx'));
    end
    w2 = [1/(mu1-x2), 1/(mu2-x2)];
    t3 = toc;
    
    figure(1), plot(N,w1(1), 'b*', 'markersize', 20); 
    hold on;
    plot(N,w1(2), 'bs', 'markersize', 20);
    plot(N,w(1), 'r*', 'markersize', 20); 
    plot(N,w(2), 'rs', 'markersize', 20);
    plot(N,w2(1), 'g*', 'markersize', 20); 
    plot(N,w2(2), 'gs', 'markersize', 20);

    figure(2), plot(N,x1, 'b*', 'markersize', 20); 
    hold on;
    plot(N,x(1), 'r*', 'markersize', 20);
    plot(N,x2, 'g*', 'markersize', 20);
    %    figure(3), plot(N,t1, 'b*', 'markersize', 20); 
    %    hold on;
    %    plot(N,t2, 'r*', 'markersize', 20);
    %    plot(N,t3, 'g*', 'markersize', 20);
    figure(4), plot(N, p, 'b*', 'markersize', 20);
    hold on;
end

figure(1), title('Waiting Times');
figure(2), title('Throughput');
%figure(3), title('Computational Time');
figure(4), title('Feedback probability');