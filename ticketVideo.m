% This script generates a video of debris pickup activity over the
% duration of a mission
% 
% 30 August 2011
% J.Brooks

% Need to run debrisanalysis.m to load a dataset
startTime = 40671; %min([towerTime; loadTime]);
endTime = max([towerTime; loadTime]);

% AVI file setup
aviobj = avifile('alDataCap.avi');
%aviobj.Quality = 100;
aviobj.Fps = 10;

% Parameters
INC = 5/60/24;                         % 5-min increments of time
FADE = 30/60/24;                       % fade time (30 min)
MAX_SIZE = 60;                          % initial size 
MIN_SIZE = 1;                           % end size
DELTA_SIZE = (MAX_SIZE-MIN_SIZE)*(INC/FADE);

active.lat = [];                       % active dataset
active.lon = [];
active.size = [];
active.cap = [];

% Create figure with state outline
figure(1);
%states = shaperead('usastatehi', 'UseGeoCoords', true);
%AL = states(1);
AL = shaperead('/home/james/Documents/Research/Debris Task Assignment/GIS/ALcounties/tl_2010_01_cousub10.shp', 'UseGeoCoords', true);
colormap(jet(50));
caxis([10 120]);
caxis('Manual');
colorbar('East');
%geoshow(AL); hold on; 
set(gcf, 'Position', [1 1 1280 721]);

for i = 1:length(AL)
    if sum(AL(i).BoundingBox(:,2) > 32) > 1
        plot(AL(i).Lon, AL(i).Lat, 'b'); hold on;    
    end
end

axis equal;

for i = startTime: INC: endTime

    if ~isempty(active.size)
        active.size = active.size - DELTA_SIZE;
        keepTickets = find(active.size >= MIN_SIZE);
        active.lat = active.lat(keepTickets);
        active.lon = active.lon(keepTickets);
        active.size = active.size(keepTickets);
        active.cap = active.cap(keepTickets);
        delete(tmp);
    end
        
    newTickets = find(loadTime > i & loadTime < i+INC & lat~=0 & lon~=0);
    
    if ~isempty(newTickets)
        newCap = zeros(length(newTickets), 1);
        for j = 1:length(newTickets)
            newCap(j) = capacity(find(trucks == truckId(newTickets(j))));
        end
        active.lat = [active.lat; lat(newTickets)];
        active.lon = [active.lon; lon(newTickets)];
        active.size = [active.size; repmat(MAX_SIZE, ...
                                           length(newTickets),1)];
        active.cap = [active.cap; newCap];
    end

    tmp = scatter(active.lon, active.lat, active.size, active.cap, 'filled');
    title(sprintf('Days: %2.2f', i - startTime));
    caxis([10 120]);
    colorbar('East');    
    aviobj = addframe(aviobj,getframe());
end
aviobj = close(aviobj);