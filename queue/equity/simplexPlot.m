all close;
revGray = flipud(colormap('gray'));

figure, surf(w1,w2,w3,throughput,'EdgeColor','none');
xlabel('w_1'); ylabel('w_2'); zlabel('w_3');
title('Throughput')
colorbar;
view([0.5 0.5 0.5])
hold on, plot3([1 0 0 1], [0 1 0 0], [0 0 1 0], 'k');
colormap(gcf, gray);
caxis([0, 2.2]);

% --------------------
% PLOT SIMPLEX
% --------------------

figure, surf(w1,w2,w3,inequity/sqrt(6),'EdgeColor','none');
xlabel('w_1'); ylabel('w_2'); zlabel('w_3');
title('Inequity - Flows');
colorbar;
view([0.5 0.5 0.5])
hold on, plot3([1 0 0 1], [0 1 0 0], [0 0 1 0], 'k');
colormap(gcf, revGray);
caxis([0, 1]);

figure, surf(w1,w2,w3,inequity2/sqrt(6),'EdgeColor','none');
xlabel('w_1'); ylabel('w_2'); zlabel('w_3');
title('Inequity - Wait Times');
colorbar;
view([0.5 0.5 0.5])
hold on, plot3([1 0 0 1], [0 1 0 0], [0 0 1 0], 'k');
colormap(gcf, revGray);
caxis([0, 1]);

figure, surf(w1,w2,w3,inequity3/sqrt(6),'EdgeColor','none');
xlabel('w_1'); ylabel('w_2'); zlabel('w_3');
title('Inequity - Lengths');
colorbar;
view([0.5 0.5 0.5])
hold on, plot3([1 0 0 1], [0 1 0 0], [0 0 1 0], 'k');
colormap(gcf, revGray);
caxis([0, 1]);

figure, surf(w1,w2,w3,inequity4/sqrt(6),'EdgeColor','none');
xlabel('w_1'); ylabel('w_2'); zlabel('w_3');
title('Inequity - Utilization');
colorbar;
view([0.5 0.5 0.5])
hold on, plot3([1 0 0 1], [0 1 0 0], [0 0 1 0], 'k');
colormap(gcf, revGray);
caxis([0, 1]);

% --------------------
% PLOT FRONTIER
% --------------------

figure, plot(reshape(throughput,1,N*N), reshape(inequity/sqrt(6),1,N*N), 'k.');
xlabel('Throughput');
ylabel('Inequity');
title('Efficient Frontier - Flows');

figure, plot(reshape(throughput,1,N*N), reshape(inequity2/sqrt(6),1,N*N), 'k.');
xlabel('Throughput');
ylabel('Inequity');
title('Efficient Frontier - Wait Times');

figure, plot(reshape(throughput,1,N*N), reshape(inequity3/sqrt(6),1,N*N), 'k.');
xlabel('Throughput');
ylabel('Inequity');
title('Efficient Frontier - Lengths');

figure, plot(reshape(throughput,1,N*N), reshape(inequity4/sqrt(6),1,N*N), 'k.');
xlabel('Throughput');
ylabel('Inequity');
title('Efficient Frontier - Utilization');