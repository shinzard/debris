%data = csvread('cogBeh_WTC_OKC.csv');
% data = csvread('allCogBehCodes_2oct2011.csv');
data = csvread('allCogBehCodes_7nov2011_self.csv');

% count self data only
self = data(:,6);
keep = find(self == 1);

figure, hist(data(keep,5), [1:8])

docs = unique(data(keep,1)+data(keep,3)*1000);
doc = data(keep,1) + data(keep,3)*1000;
dis = data(keep,3);
sent = data(keep,2);
code = data(keep,5);
figure, plot(code, 'b.')
type = data(keep,4);

numSent = zeros(1,length(docs));

for i = 1:length(docs)
    idx = find(doc == docs(i));
    numSent(i) = max(sent(idx));
end

figure, hist(numSent, 20)
title('Histogram of Number of Sentences Per File')
xlabel('Number of Sentences')
ylabel('Document Count')
text(100, 4, sprintf('Mean Sentences: %2.2f', mean(numSent)))
text(100, 4, sprintf('Total Sentences: %2.2f', sum(numSent)))

axis([0 250 0 8])

for i = 1:length(docs)
%     cognitive count excludes intervention
idxC = find(doc == docs(i) & type == 1 & code < 7 & code > 1);
idxB = find(doc == docs(i) & type == 2 & code < 7);
numC(i) = length(idxC);
numB(i) = length(idxB);
end

figure, hist(numB, 20)
xlabel('Number of Behavioral Events')
ylabel('Document Count')
title('Histogram of Behavioral Codes per Document')
text(20, 10, sprintf('Mean BEH: %2.2f', mean(numB)))
text(20, 9, sprintf('Total BEH: %d', sum(numB)))

figure, hist(numC, 20)
title('Histogram of Cognitive Codes per Document')
xlabel('Number of Cognitive Events')
ylabel('Document Count')
text(150, 15, sprintf('Mean COG: %2.2f', mean(numC)))
text(150, 13, sprintf('Total COG: %d', sum(numC)))
