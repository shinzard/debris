clear all;
close all;

load DEBRIS_AL_COMPLETE

projs = unique(project);

for i = 1:length(projs)
    figure, plot(lon,lat, 'b.'), hold on;
    idx = find(project == projs(i));
    plot(lon(idx),lat(idx), 'r.');
    title(sprintf('Project %2.1f',projs(i)));
    axis([-88.5,-85,32,35.3])
end

placefigures;