function struc = groupId(truckId, QC, loadTime, subcont, haulMi, ...
                         MIN_CORR, PLOTS)
% GROUPID - This function calculates the most probable group assignments for a debris removal mission
%   truckId
%   QC

% 10 September 2011
% J.Brooks

% LAST WORKING ON GETTING GROUP ID's TO SPAN ADJACENT DAYS....NOT COMPLETE
    
    if nargin < 6
        PLOTS = 0;
        MIN_CORR = 0.9;
    end
    
    day = [min(floor(loadTime)):max(floor(loadTime))];
    QCs = unique(QC);
    trucks = unique(truckId);
    
    % Day stats
    struc.numCrews = zeros(1,length(day));
    struc.crewSize = [];
    struc.percSpan = zeros(1,length(day));
    struc.avgSize = zeros(1,length(day));
    struc.sdSize = zeros(1,length(day));
    struc.minSize = zeros(1,length(day));
    struc.maxSize = zeros(1,length(day));

    % Group ID for each ticket
    struc.GroupID = zeros(1,length(truckId));
    groupId = 1;

    
    for i = 5:length(day)
        k = 1;
        disp(sprintf('Day %d', day(i)));
        
        % Clear structure
        clear groupStruc;
        groupStruc.Trucks = [];
        groupStruc.QC = [];
        groupStruc.Sub = [];
        groupStruc.NumTickets = [];
        groupStruc.HaulAvg = [];
        groupStruc.HaulSD = [];
        groupStruc.ID = [];
        if exist('yGroupID')
            lastYgroupID = yGroupID;
        end
        yGroupID = zeros(1,length(QCs));
        
    
        % Find day data
        idx = find(floor(loadTime) == day(i) & QC ~= 0 & ...
                   subcont ~= 0 & haulMi < 56);
        
        if ~isempty(idx)
        
            % Store yesterday
            if exist('table')
                yesterday = table;
            end

            % Find today's correlation between QCs
            table = contingencyTable([truckId(idx), QC(idx)], trucks, ...
                                     QCs, 't', 'q', 10, 0);

            if exist('yesterday') 
                
                % Find correlation with yesterday's data
                yRho = corr(yesterday>0, table>0);

                if PLOTS
                    figure, imagesc(yRho);
                    title(sprintf(['Yesterday Correlation QC/Truck ' ...
                                   'Assignment; Day: %d'], i));
                    colormap(bone(5))
                    caxis([0 1])
                end
            else
                yRho = zeros(size(table));
            end
            
            %truL = unique(truckId(idx));
            %qcL = unique(QC(idx));
            
            rho = corr(table > 0);
            highC = rho > MIN_CORR;
            
            remaining = ones(1, size(rho,1));
            
            for j = 1:size(rho,1)
                % Find highly correlated QC groups that haven't
                % been combined already
                qcTmp = find( (highC(j,:) & remaining) > 0);
                
                if ~isempty(qcTmp)
                    [t, q] = ind2sub(size(table(:,qcTmp)), ...
                                     find(table(:,qcTmp) > 0));
                    
                    groupStruc(k).Trucks = trucks(t);
                    groupStruc(k).QC = QCs(qcTmp);

                    % Assign unique group ID
                    yCor = find(yRho(j,:) > MIN_CORR);
                    if isempty(yCor)
                        groupStruc(k).ID = groupId;
                        yGroupID(j) = groupId;
                        groupId = groupId + 1;
                    else
                        groupStruc(k).ID = lastYgroupID(yCor);
                    end

                    groupTickets = [];
                    for l = 1:length(t)
                        for m = 1:length(qcTmp)
                            groupTickets = [groupTickets, 
                                  find(truckId(idx) == trucks(t(l)) & ...
                                  QC(idx) == QCs(qcTmp(m)))];
                        end
                    end
                    groupStruc(k).Sub = ...
                        unique(subcont(idx(groupTickets)));
                    groupStruc(k).HaulAvg = ...
                        mean(haulMi(idx(groupTickets)));
                    groupStruc(k).HaulSD = ...
                        std(haulMi(idx(groupTickets)));
                    groupStruc(k).NumTickets = ...
                        length(groupTickets);
                    remaining(qcTmp) = 0;
                    k = k + 1;
                end
            end
            
            for j = 1:length(groupStruc)
                groupStruc(j).Size = length(groupStruc(j).Trucks);
                groupStruc(j).Span = length(groupStruc(j).Sub);
            end
            
            struc.crewSize = [struc.crewSize, collapse(groupStruc, 'Size')];
            % capture statistics
            struc.numCrews(i) = length(groupStruc);
            disp(sprintf('Number of crews: %d', struc.numCrews(i)));
            struc.percSpan(i) = length(find(collapse(groupStruc, ...
                                 'Span') > 1))/length(groupStruc); 
            struc.avgSize(i) = mean(collapse(groupStruc, 'Size'));
            struc.sdSize(i)  = std(collapse(groupStruc, 'Size'));
            struc.minSize(i)  = min(collapse(groupStruc, 'Size'));
            struc.maxSize(i)  = max(collapse(groupStruc, 'Size'));
            
            if PLOTS
                figure(1), plot(collapse(groupStruc,'Size'), ...
                                collapse(groupStruc,'NumTickets')./ ...
                                collapse(groupStruc,'Size'), 'b.');
                hold on;
                xlabel('Crew Size (# trucks)');
                ylabel('Number of tickets');
                title('Crew Productivity by Size');
                
                figure(2), plot(collapse(groupStruc,'Size'), ...
                                collapse(groupStruc,'HaulAvg'), 'b.', ...
                                collapse(groupStruc,'Size'), ...
                                collapse(groupStruc,'HaulSD'), 'r.');
                hold on;
                xlabel('Crew Size (# trucks)');
                ylabel('Number of tickets');
                title('Crew Productivity by Size');
            end
        end        
    end    
    
    struc.Last = groupStruc;
    
    figure, plot(struc.numCrews), title('Number of Crews')
    xlabel('Time (days)')
    ylabel('Crews')
    
    figure, plot(struc.percSpan)
    title('Percentage of Crews spanning 2 subcontractors');
    xlabel('Time (days)')
    ylabel('Percent')
    
    figure, plot(struc.avgSize), title('Average Crew Size')
    xlabel('Time (days)')
    ylabel('Crew size (# trucks)')