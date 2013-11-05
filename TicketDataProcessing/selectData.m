function teams = selectData(tickets,duration)
    
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