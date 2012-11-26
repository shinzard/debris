close all;
clear all;

mu1 = 8;                                % *
mu2 = 7;                                % s
mu3 = 20;                               % ^

P = [0,0,1; 0,0,1 ; 0.5,0.5,0];
V = qnvisits(P);
S = [1/mu1, 1/mu2, 1/mu3];

% For new approx method
l = 0.1;

for N = 1:30;
    % New Approx Method
    tic;                                % can we improve with lower
                                        % lambda (higher p)??
    
    % These terms occur often...    
    l1 = l*(N+1);
    l2 = l^2*(N+2);
    l3 = l^3*(N+3);
    prod = N*mu1*mu2*(mu3+l);
    
    pRoot = roots([-(l2 + 2*l1*(mu1+mu2)+4*prod/(mu3+l))*(mu3+l), ...
                   -(l3 + l2*(2*mu1+2*mu2-(mu3+l)) + 4*l1*(mu1*(mu2-(mu3+l))-mu2*(mu3+l))-12*prod), ...
                   2*(l2*(mu1+mu2)+l1*(mu1*(4*mu2-(mu3+l))-mu2*(mu3+l)) - 6*prod),...
                   -4*(l1-N*(mu3+l))*mu1*mu2]);
    p = min(pRoot);
    w1 = [ 2*(1-p)/(2*(1-p)*mu1 - p*l), 2*(1-p)/(2*(1-p)*mu2 - p*l), (1-p)/((1-p)*(mu3+l) - l)];
    x1 = l/(1-p);
    t1 = toc;

    % check
    %    l1 = l/((1-p)*mu1 - l); 
    %    l2 = p*l/((1-p)*mu2 - p*l);
    %    l1+l2

    % MVA
    tic;
    [u,w,q,x] = qnclosed(N,S,V);
    t2 = toc;
    
    
    figure(1), plot(N,w1(1), 'b*', 'markersize', 20); 
    hold on;
    plot(N,w1(2), 'bs', 'markersize', 20);
    plot(N,w1(3), 'b^', 'markersize', 20);
    plot(N,w(1), 'r*', 'markersize', 20); 
    plot(N,w(2), 'rs', 'markersize', 20);
    plot(N,w(3), 'r^', 'markersize', 20);

    figure(2), plot(N,x1, 'b*', 'markersize', 20); 
    hold on;
    plot(N,x(3), 'r*', 'markersize', 20);

    figure(4), plot(N, p, 'b*', 'markersize', 20);
    hold on;
end

figure(1), title('Waiting Times');
figure(2), title('Throughput');
%figure(3), title('Computational Time');
figure(4), title('Feedback probability');