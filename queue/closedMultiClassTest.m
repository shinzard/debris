## closedTest

## Author: james <james@brooks>
## Created: 2012-03-20

function throughput = closedMultiClassTest(p2)
    
close all;

P = zeros(2,3,2,3); %[0, 0.3, 0.7; 1, 0, 0; 1, 0, 0];
P(1,1,1,3) = 1;
P(1,3,1,1) = 1;

P(2,1,2,2) = 1;
P(2,2,2,1) = 1;

S = [0.1, 0.2, 0.2;0.1, 0.2, 0.2];
V = qnvisits(P);

%U = zeros(3,2,100);
%R = zeros(3,2,100);
%Q = zeros(3,2,100);
%X = zeros(3,2,100);
throughput = zeros(10,10);
length1 = zeros(10,10);
length2 = zeros(10,10);
length3 = zeros(10,10);

for N1 = 1:10
    for N2 = 1:10
        N = [N1, N2]; %[round(total*p2), total - round(total*p2)]
        [u,r,q,x] = qnclosedmultimva(N,S,V);
        throughput(N1,N2) = x(1,1)/V(1,1) + x(2,1)/V(2,1);
        length1(N1,N2) = sum(q(:,1));
        length2(N1,N2) = sum(q(:,2));
        length3(N1,N2) = sum(q(:,3));
    end
end

%figure, plot([1:100], U); title('Utilization');
%figure, plot([1:100], R); title('Response Time');
%figure, plot([1:100], Q); title('Avg Num Customers');
%figure, plot([1:100], X); title('Throughput');

figure, mesh([1:10],[1:10], throughput);
xlabel('Customers Assigned to Server 2');
ylabel('Customers Assigned to Server 3');
title('System Throughput')

figure, mesh([1:10],[1:10], length1);
xlabel('Customers Assigned to Server 2');
ylabel('Customers Assigned to Server 3');
title('Mean Queue Length: Server 1')

figure, mesh([1:10],[1:10], length2);
xlabel('Customers Assigned to Server 2');
ylabel('Customers Assigned to Server 3');
title('Mean Queue Length: Server 2')

figure, mesh([1:10],[1:10], length3);
xlabel('Customers Assigned to Server 2');
ylabel('Customers Assigned to Server 3');
title('Mean Queue Length: Server 3')