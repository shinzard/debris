close all;
clear all;
more off;

% DoE
mu1 = 10;
mu2V = 20.8; %[1:2:25];
mu3V = 45; %[40:5:80];
NV = 10; %[6:6:24];

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
poa = zeros(NUM_TESTS,1);

test = 1;
plotIdx = 1;

for mu3 = mu3V
    disp(sprintf('Running mu3 = %d', mu3));
    for mu2 = mu2V
        disp(sprintf('Running mu2 = %d', mu2));
        for N = NV
            opt = 0;
            mu3r(test) = mu3; mu2r(test) = mu2; nr(test) = N;
            S = [1/mu3, 1/mu1, 1/mu2];

            % Rational equilibrium
            [ut,wt,qt,xt] = queueSim(N,S,-1,3000,500)
            
            % conjectured optimal
            [ut1,wt1,qt1,xt1] = queueSim(N,S,-2,3000,500)
            
            poa(test) = xt(1)/xt1(1);
            disp(sprintf('Price of Anarchy: %2.2f', poa(test)));
            test = test + 1;
        end
    end
end

figure, plot(test, poa, 'b.');
title('Price of Anarchy');