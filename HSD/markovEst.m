% This function estimates the markov chain underlying 
% cognition events
% 
% 2 Oct 2011
% J.Brooks

beh = find(type == 2 & code < 7);

% Initialize sequence matrix
seq = zeros(length(beh), 7);

for i = 1:length(beh)
    idx = find(doc == doc(i) & sent == sent(i) & type == 1 & code < 7);
    seq(i,1:length(idx)) = code(idx);
end

na3 = find(seq(:,3) == 0);
e3 = find(seq(:,3) == 1);
h3 = find(seq(:,3) == 2);
o3 = find(seq(:,3) == 3);
g3 = find(seq(:,3) == 4);

na2 = find(seq(:,2) == 0);
e2 = find(seq(:,2) == 1);
h2 = find(seq(:,2) == 2);
o2 = find(seq(:,2) == 3);
g2 = find(seq(:,2) == 4);

na1 = find(seq(:,1) == 0);
e1 = find(seq(:,1) == 1);
h1 = find(seq(:,1) == 2);
o1 = find(seq(:,1) == 3);
g1 = find(seq(:,1) == 4);

prob = zeros(21, 4);

n =  hist(seq(:, 3), [0:4])
prob(1:5,1) = n/sum(n);

%% First stage
n =  hist(seq(na3, 2), [0:4]);
prob(1:5,2) = n/sum(n);

n =  hist(seq(e3, 2), [1:4]);
prob(6:9,2) = n/sum(n);

n =  hist(seq(h3, 2), [1:4]);
prob(10:13,2) = n/sum(n);

n =  hist(seq(o3, 2), [1:4]);
prob(14:17,2) = n/sum(n);

n =  hist(seq(g3, 2), [1:4]);
prob(18:21,2) = n/sum(n);


%% Second stage
n =  hist(seq(na2, 1), [0:4]);
prob(1:5,3) = n/sum(n);

n =  hist(seq(e2, 1), [1:4]);
prob(6:9,3) = n/sum(n);

n =  hist(seq(h2, 1), [1:4]);
prob(10:13,3) = n/sum(n);

n =  hist(seq(o2, 1), [1:4]);
prob(14:17,3) = n/sum(n);

n =  hist(seq(g2, 1), [1:4]);
prob(18:21,3) = n/sum(n);

%% Final stage
n = hist(code(beh(na3)), [5:6]);
prob(1:2, 4) = n/sum(n);

n = hist(code(beh(e3)), [5:6]);
prob(6:7, 4) = n/sum(n);

n = hist(code(beh(h3)), [5:6]);
prob(10:11, 4) = n/sum(n);

n = hist(code(beh(o3)), [5:6]);
prob(14:15, 4) = n/sum(n);

n = hist(code(beh(g3)), [5:6]);
prob(18:19, 4) = n/sum(n);