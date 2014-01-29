C = 1;
R = 5;
mu = [4 6];

maxT = R - C*(1/4+1/6);
Tvec = linspace(0,maxT,1000);
%Pvec = linspace(0.01,0.99,20);
% maxToll = zeros(length(Pvec),length(Tvec));
% maxSToll = zeros(length(Pvec),length(Tvec));
% S = zeros(length(Pvec),length(Tvec));
% I = zeros(length(Pvec),length(Tvec));
% thruRatio = zeros(length(Pvec),length(Tvec));
tollRate = zeros(1,length(Tvec));
ns = zeros(1,length(Tvec));
indivRate = zeros(1,length(Tvec));
thru = zeros(1,length(Tvec));

i = 1;
for T = Tvec
    %        disp(sprintf('T = %f',T));
    [tollRate(i), ns(i), indivRate(i), thru(i)] = explore(C,R,T,mu,0);
    i = i + 1;
end

close all;
disp('done!')
%mesh(Tvec, Pvec, maxToll)
plot(Tvec, tollRate, 'b');
set(gca, 'fontsize', 20);
hold on;
plot(Tvec, ns, 'k');
plot(Tvec, indivRate, 'b--');
plot(Tvec, indivRate + tollRate, 'm');
plot(Tvec, thru, 'g');
title('Toll, Individual Rate and Threshold');
xlabel('Per-Cycle Toll Levied');
legend({'Toll Rate', 'Joining Threshold', 'Indiv. Earning Rate', ...
        'Overall Earning Rate', 'System Throughput'}, 'location', 'northeast');

[val, idx] = max(tollRate);
disp(sprintf('Max Toll: %f \t Rate: %f \t Threshold: %d \t Indiv: %f', Tvec(idx), ...
             val, ns(idx), indivRate(idx)));