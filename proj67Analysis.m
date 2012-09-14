load DEBRIS_ANALSIS_AL_COMPLETE
idx = find(project == 6 | project == 7);

 
tdsr602 = find(tdsr(idx)==602);
tdsr601 = find(tdsr(idx)==601);
t602 = sort(towerTime(idx(tdsr602)));
t601 = sort(towerTime(idx(tdsr601)));

% Remove truck/trailor (18 second cutoff)
idx = find([1, diff(t601)'*24]>0.005); 
t601 = t601(idx);

%%
plotArrivals(t601, floor(t601))
%plotArrivals(t602, floor(t602))
%%
days = unique(floor(t601)');
for i = days
    idx = find(floor(t601)==i);
    if length(idx)>150
        figure, hist(diff(t601(idx)*24),100);
        title(sprintf('Day: %d',i));
    end
end
%% 
group = QC(idx)*1e6+floor(loadTime(idx));
groups = unique(group);

colors = ['brkgmyc'];

for i = 1:length(groups)
   c = colors(mod(i-1,length(colors))+1);
   tmp = find(group == groups(i));
   tdsrTmp = tdsr(idx(tmp));
   tdsrTmps = unique(tdsr(idx(tmp)));
   %disp(sprintf('Group: %d, Num TDSRs: %d', groups(i), length(unique(tdsrTmp))));
   for j = 1:length(unique(tdsrTmp));
       tmp2 = find(tdsrTmp == tdsrTmps(j));
       % one-way (loaded) travel time
       if tdsrTmps(j)==601
        figure(1), plot(haulMi(idx(tmp(tmp2))), (towerTime(idx(tmp(tmp2)))-loadTime(idx(tmp(tmp2))))*24, [c,'.']);
        hold on;
       else
        figure(2), plot(haulMi(idx(tmp(tmp2))), (towerTime(idx(tmp(tmp2)))-loadTime(idx(tmp(tmp2))))*24, [c,'.']);
        hold on;
           
       %pause(2);
       end
   end
   
   figure(1), title('601'), figure(2), title('602');
end