% Finite-Source
[u,w,q,x,e]=optimalAssignment(mu,type,0,20,3)
figure(1), 
subplot(411), plot(u,'r.'), hold on;
subplot(412), plot(w,'r.'), hold on;
subplot(413), plot(q,'r.'), hold on;
subplot(414), plot(x,'r.'), hold on;

subplot(411), title('Utilization');
subplot(412), title('Wait Times');
subplot(413), title('Length')
subplot(414), title('Flows');

% M/M/1/N
[u,w,q,x,e]=optimalAssignment(mu,type,0,20,2);
subplot(411), plot(u,'g.');
subplot(412), plot(w,'g.');
subplot(413), plot(q,'g.');
subplot(414), plot(x,'g.');

% M/M/1
[u,w,q,x,e]=optimalAssignment(mu,type,0,20,1);
subplot(411), plot(u,'b.');
subplot(412), plot(w,'b.');
subplot(413), plot(q,'b.');
subplot(414), plot(x,'b.');

legend({'finite-source', 'mm1k', 'mm1'});