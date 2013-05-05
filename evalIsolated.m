clear all;
close all;

load DEBRIS_AL_COMPLETE

disp(sprintf('>15\t # subs\t # tickets'));
idx = find(project == 6 | project == 7);
disp(sprintf('%2.2f\t %d\t %d', sum(haulMi(idx)>15)/length(idx), ...
    length(unique(subcont(idx))), length(idx)));
x = match(unique(subcont(idx)),subcont(idx));
figure, subplot(211), hist(x,max(x));
x = match(unique(tdsr(idx)),tdsr(idx));
subplot(212), hist(x),max(x);


idx = find(project == 8 | project == 9 | project == 10);
disp(sprintf('%2.2f\t %d\t %d', sum(haulMi(idx)>15)/length(idx), ...
    length(unique(subcont(idx))), length(idx)));
x = match(unique(subcont(idx)),subcont(idx));
figure, subplot(211), hist(x,max(x));
x = match(unique(tdsr(idx)),tdsr(idx));
subplot(212), hist(x),max(x);


idx = find(project == 15);
disp(sprintf('%2.2f\t %d\t %d', sum(haulMi(idx)>15)/length(idx), ...
    length(unique(subcont(idx))), length(idx)));
x = match(unique(subcont(idx)),subcont(idx));
figure, subplot(211), hist(x,max(x));
x = match(unique(tdsr(idx)),tdsr(idx));
subplot(212), hist(x),max(x);


idx = find(project == 20 | project == 22);
disp(sprintf('%2.2f\t %d\t %d', sum(haulMi(idx)>15)/length(idx), ...
    length(unique(subcont(idx))), length(idx)));
x = match(unique(subcont(idx)),subcont(idx));
figure, subplot(211), hist(x,max(x));
x = match(unique(tdsr(idx)),tdsr(idx));
subplot(212), hist(x),max(x);

idx = find(project == 20);
disp(sprintf('%2.2f\t %d\t %d', sum(haulMi(idx)>15)/length(idx), ...
    length(unique(subcont(idx))), length(idx)));
x = match(unique(subcont(idx)),subcont(idx));
figure, subplot(211), hist(x,max(x));
x = match(unique(tdsr(idx)),tdsr(idx));
subplot(212), hist(x),max(x);

idx = find(project == 22);
disp(sprintf('%2.2f\t %d\t %d', sum(haulMi(idx)>15)/length(idx), ...
    length(unique(subcont(idx))), length(idx)));
x = match(unique(subcont(idx)),subcont(idx));
figure, subplot(211), hist(x,max(x));
x = match(unique(tdsr(idx)),tdsr(idx));
subplot(212), hist(x),max(x);


idx = find(project == 24 | project == 30 | project == 31);
disp(sprintf('%2.2f\t %d\t %d', sum(haulMi(idx)>15)/length(idx), ...
    length(unique(subcont(idx))), length(idx)));
x = match(unique(subcont(idx)),subcont(idx));
figure, subplot(211), hist(x,max(x));
x = match(unique(tdsr(idx)),tdsr(idx));
subplot(212), hist(x),max(x);

idx = find(project == 24);
disp(sprintf('%2.2f\t %d\t %d', sum(haulMi(idx)>15)/length(idx), ...
    length(unique(subcont(idx))), length(idx)));
x = match(unique(subcont(idx)),subcont(idx));
figure, subplot(211), hist(x,max(x));
x = match(unique(tdsr(idx)),tdsr(idx));
subplot(212), hist(x),max(x);

placefigures;