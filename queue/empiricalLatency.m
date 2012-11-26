close all;
clear all;
more off;

mu1 = 10;
mu2V = [1:2:25];
mu3V = [5:5:40];
NV = [6:3:24];

% System configuration
P = zeros(2,3,2,3);
P(1,1,1,3) = 1;
P(1,3,1,1) = 1;
P(2,2,2,3) = 1;
P(2,3,2,2) = 1;
V = qnvisits(P);

% Results vectors
NUM_TESTS = length(mu2V)*length(mu3V)*length(NV);
lambda = zeros(NUM_TESTS*3,1);
mu = zeros(NUM_TESTS*3,1);
L = zeros(NUM_TESTS*3,1);

test = 1;

for mu3 = mu3V
    disp(sprintf('Running mu3 = %d', mu3));
    for mu2 = mu2V
        disp(sprintf('Running mu2 = %d', mu2));
        for N = NV
            opt = 0;
            mu3r(test) = mu3; mu2r(test) = mu2; nr(test) = N;
            S2 = repmat([1/mu1, 1/mu2, 1/mu3],2,1);

            % look at all possible pure assignments
            for N1 = [1:N-1]
                [u,w,q,x] = qnclosed([N1,N-N1],S2,V);
                lambda(3*(test-1)+1:3*test) = [x(1,1), x(2,2), sum(x(:,3))];
                mu(3*(test-1)+1:3*test) = [mu1, mu2, mu3];
                L(3*(test-1)+1:3*test) = [q(1,1), q(2,2), sum(q(:,3))];
            end

            test = test + 1;
        end
    end
end

