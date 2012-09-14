function struc = disagreement(data, text1, text2)
% function struc = disagreement(data, text1, text2)
% 
% This function generates levels of disagreement in data betweentwo
% categorical variables contained in columns 1 and 2 of the data
% input. Text1(2) are optional axis labels input. 
% 
% 2 September 2011
% J.Brooks
    
    if nargin < 3
        text1 = 'cat1';
        text2 = 'cat2';
    end
    
    try
        if ~isvector(data)
            data = reshape(data, max(size(data)),2);
        else
            data = reshape(data, 1, 2);
        end
    catch
        error('Only two categories can be included');
    end
    
    cat1 = unique(data(:,1));
    cat2 = unique(data(:,2));

    % Allocate structure (how?)
    
    for i = 1:length(cat1)
        idx = find(data(:,1) == cat1(i));
        cat2Idx = match(cat2, data(idx,2));
        
        struc(i).Counts = hist(cat2Idx, [1:length(cat2)]);
        [maxCnt, ID] = max(struc(i).Counts);
        struc(i).MajorityID = ID;
        struc(i).MajorityPerc = maxCnt/sum(struc(i).Counts);
    end
    