function [teams,sub] = selectData(tickets,duration)
    
    QCday = tickets.QC + 1e6*floor(tickets.loadTime);
    
    if nargin == 1
        teams = unique(QCday);
    end
    
    if nargin == 2
        idx = [];
        for i = 1:length(duration)
            idx = [idx; find(floor(tickets.loadTime)==duration(i))];
        end
        teams = unique(QCday(idx));
    end
    
    disp(sprintf('Num unique QCdays: %d', length(teams)));
    
    % Get associated subcontractor 
    sub = zeros(1,length(teams));
    subcont = match(unique(tickets.subcont), tickets.subcont,999);
    for i = 1:length(teams)
        idx = find(QCday == teams(i));
        subs = unique(subcont(idx));
        if length(subs) == 1
            sub(i) = subs;
        else
            sub(i) = NaN;
        end
    end