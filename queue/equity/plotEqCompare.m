

figure, plot(Imu, waitEq./opt, 'ks', 'MarkerSize', 10);
hold on;
semilogx(Imu, lengthEq./opt, 'k*');
semilogx(Imu, flowEq./opt, 'k^', 'MarkerSize', 8);
title('Comparison of Effectiveness By Equity Measure');
xlabel('Inequity Measure of Service Rates');
ylabel('Relative Total Throughput');

legend({'Wait', 'Length', 'Flow'})

figure, plot(optTOconvex, ItoConvex, 'b.');
hold on;
plot(flowEq, zeros(1,length(flowEq)), 'rs');
title('Eff-Eq Tradeoff Midpoint');
xlabel('System Throughput');
ylabel('Flow Inequity');
