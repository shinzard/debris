function plotDayData(qcSpan, loadTime, truckId, QC, project, lat, lon, dayData, day, region)
% PLOTDAYDATA - Display daily data
%   Generate one plot for each day in the input vector. If this
%   vector is omitted, the entire mission is shown.
%
% 6 Sept 2011
% J.Brooks

    close all;
    
    if nargin < 9
        day = [floor(min(loadTime)): floor(max(loadTime))];
        region = 'all';
    end
               
    % Collapse QC data
    qclatP = collapseV(qcSpan, 'latP');
    qclonP = collapseV(qcSpan, 'lonP');
    qclonSpan = collapseV(qcSpan, 'lonSpan');
    qclatSpan = collapseV(qcSpan, 'latSpan');
    [qcDay, qcidx] = collapseV(qcSpan, 'day');
    QCs = unique(QC);
    
    % Generate plot for each day
    for j = day
        % Find day data
        try
            idx = find(floor(loadTime) == j & lat~=0 & lon~=0 & project ...
                       == region);
        catch
            idx = find(floor(loadTime) == j & lat~=0 & lon~=0);
        end
        activeQc = match(QCs, unique(QC(idx)));
        tmp = find(~isnan(match(activeQc, qcidx)));
        qcI = find(qcDay(tmp) == j);
        
        
        % Find correlation between QCs
        table = contingencyTable([truckId(idx), QC(idx)]);
        rho = corr(table>0);
        
        figure, plot(lat(idx), lon(idx), 'b.')
        for i = 1:length(qcI)
            obj = rectangle('Position', [qclatP(tmp(qcI(i))), ...
                                qclonP(tmp(qcI(i))), ...
                                max(0.0002, qclatSpan(tmp(qcI(i)))), ...
                                max(0.0002, qclonSpan(tmp(qcI(i))))], ...
                            'Curvature', [0.05, 0.05]);
            if length(find(rho(i,:) > 0.8)) > 1
                set(obj, 'EdgeColor', [0.8 0 0]);
            end
        end
        % [cluster, centers, sum, distance] = kmeans([lat(idx), lon(idx)], 159);
        % plot(centers(:,1), centers(:,2), 'r.', 'MarkerSize', 6)
        
        title(sprintf('Day: %d', j));
        
    end
