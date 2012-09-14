## closedTest

## Author: james <james@brooks>
## Created: 2012-03-20

close all;
clear all;

P = [0, 0.3, 0.7; 1, 0, 0; 1, 0, 0];
S = [10, 30, 30];
V = qnvisits(P);

total = 10;

i = 1;
U = zeros(3,total);
R = zeros(3,total);
Q = zeros(3,total);
X = zeros(3,total);

for N = 1:total
    disp(sprintf('N: %d', N));
    [U(:,i), R(:,i), Q(:,i), X(:,i)] = qnclosed(N,S,V);
    i = i + 1;
end

figure, plot([1:total], U); title('Utilization');
figure, plot([1:total], R); title('Response Time');
figure, plot([1:total], Q); title('Avg Num Customers');
figure, plot([1:total], X); title('Throughput');