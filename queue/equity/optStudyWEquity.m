close all;
clear all;
more off;

mu = 4;
mu1V = [0.1:0.5:4];
mu2V = [0.1:0.5:4];
mu3V = [0.1:0.5:4];
N = 5; %[6:3:24];

% Results vectors
NUM_TESTS = length(mu1V)*length(mu2V)*length(mu3V)*length(N);
waitEq = zeros(1,NUM_TESTS);
lengthEq = zeros(1,NUM_TESTS);
flowEq = zeros(1,NUM_TESTS);
opt = zeros(1,NUM_TESTS);
Imu = zeros(1,NUM_TESTS);

test = 1;

for mu1 = mu1V
    for mu2 = mu2V
        for mu3 = mu3V
          mus = [mu, mu1, mu2, mu3]
          type = [1, 2, 2, 2];

          Imu(test) = 1/3*sqrt(sum((mus(2:end)./mean(mus(2:end))-1).^2));
          
          % Optimal
          [u,w,q,x,exit]=optimalAssignment(mus,type,0,N,1);
          opt(test) = sum(x);
          
          % Wait-Equity
          waitEq(test) = equityLinesearch(mus, 0, N);
          
          % Length-Equity
          lengthEq(test) = equityLinesearch(mus, mus(2:end)/sum(mus(2:end)), N);
          
          % Flow-Equity
          flowEq(test) = equityLinesearch(mus, repmat(1/3, 1, 3), N);
          
          if waitEq(test)/opt(test) > 1.01
              mus
              error('what')
          end

          test = test + 1;
          
        end
    end
end

save(['equityLineSearch-',num2str(N),'_', datestr(now)]);

plotEqCompare;