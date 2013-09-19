close all;
clear all;

load FI_simData_26Feb2013

% ------------------------------
% SAME SERVICE RATES (A)
% uniform vs. random 
% ------------------------------
boxplot(a6.unif, a6.n, 'colors', 'b', 'labels', {'','','','','','','',''}, 'sym', 'b+');
hold on;
axis manual;
boxplot(a6.rand, a6.n, 'colors', 'r');
axis([0 9 6 20]);
xlabel('Number of Trucks');
ylabel('Mean Throughput Rate');
title('System Throughput');

% ------------------------------
% SAME SERVICE RATES
% uniform: A vs. B
% ------------------------------
figure, boxplot(a6.unif, a6.n, 'colors', 'b', 'labels', {'','','','','','','',''});
hold on;
axis manual;
boxplot(b6.unif, b6.n, 'colors', 'r');
axis([0 9 6 20]);
xlabel('Number of Trucks');
ylabel('Mean Throughput Rate');
title('System Throughput');

% ------------------------------
% DIFFERENT SERVICE RATES (A)
% uniform vs. random vs. optimal
% ------------------------------
figure, boxplot(a6.unifV, a6.n, 'colors', 'b', 'labels', {'','','', ...
                    '','','','',''}, 'sym', 'b+');
hold on;
axis manual;
boxplot(a6.randV, a6.n, 'colors', 'r', 'labels', {'','','', ...
                    '','','','',''});
boxplot(a6.optimalV, a6.n, 'colors', 'g', 'sym', 'g+');
axis([0 9 6 20]);
xlabel('Number of Trucks');
ylabel('Mean Throughput Rate');
title('System Throughput');

% ------------------------------
% DIFFERENT SERVICE RATES ()
% uniform vs. random vs. optimal
% ------------------------------
figure, boxplot(b6.unifV, b6.n, 'colors', 'b', 'labels', {'','','', ...
                    '','','','',''}, 'sym', 'b+');
hold on;
axis manual;
boxplot(b6.randV, b6.n, 'colors', 'r', 'labels', {'','','', ...
                    '','','','',''});
boxplot(b6.optimalV, b6.n, 'colors', 'g', 'sym', 'g+');
axis([0 9 6 20]);
xlabel('Number of Trucks');
ylabel('Mean Throughput Rate');
title('System Throughput');

%unifAmeans = grpstats(a6.unif, a6.n);
%randAmeans = grpstats(a6.rand, a6.n);
%figure,
%plot([6:6:48], unifAmeans, 'b');
%hold on;
%plot([6:6:48], randAmeans, 'r');
