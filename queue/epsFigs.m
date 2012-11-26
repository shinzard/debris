figure(2), plot(x(:,3)), title('System Throughput');
figure(3), plot(w), title('Wwait Times');
figure(4), plot(q), title('Queue Lengths');
figure(6), plot(u), title('Utilizations');

figure(1), ylabel('Percent of N to Q1'), xlabel('\mu_1/\mu_2');
print(1,'assignment.eps')

figure(2), xlabel('\mu_1/\mu_2');
 print(2,'sysThroughput.eps')

 if (size(w,2) == 4)
 figure(3), xlabel('Test Index'), legend({'Q1', 'Q2', 'Q3', 'Central'});
print(3,'waitTime.eps')

figure(5), xlabel('Test Index'), legend({'Q1', 'Q2', 'Q3', 'Central'}, 'location', 'northwest');
print(5,'throughput.eps')

figure(4), xlabel('Test Index'), legend({'Q1', 'Q2', 'Q3', 'Central'}, 'location', 'northeast');
print(4,'lengths.eps')

figure(6), xlabel('Test Index'), legend({'Q1', 'Q2', 'Q3', 'Central'}, 'location', 'northeast');
print(6,'utilizations.eps')

 end
 
 
 figure(3), xlabel('Test Index'), legend({'Q1', 'Q2', 'Central'});
print(3,'waitTime.eps')

figure(5), xlabel('Test Index'), legend({'Q1', 'Q2', 'Central'}, 'location', 'northwest');
print(5,'throughput.eps')

figure(4), xlabel('Test Index'), legend({'Q1', 'Q2', 'Central'}, 'location', 'northeast');
print(4,'lengths.eps')

figure(6), xlabel('Test Index'), legend({'Q1', 'Q2', 'Central'}, 'location', 'northeast');
print(6,'utilizations.eps')
