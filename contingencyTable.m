function counts = contingencyTable(data, cat1, cat2, text1, text2, maxC, PLOTS)
% function counts = contingencyTable(data, text1, text2, maxC)
% 
% This function generates contingency table to be used for chi-squared
% tests of independence between two categorical variables contained
% in columns 1 and 2 of the data input. Text1(2) are optional axis
% labels input. maxC sets colormap maximum.
% 
% 31 August 2011
% J.Brooks
    
    if nargin < 3
        cat1 = unique(data(:,1));
        cat2 = unique(data(:,2));
        text1 = 'cat1';
        text2 = 'cat2';
        PLOTS = 0;
    end
    
    try
        if ~isvector(data)
            data = reshape(data, max(size(data)), 2);
        else
            data = reshape(data, 1, 2);
        end
    catch
        error('Only two categories can be included');
    end
    
    % Allocate table
    counts = zeros(length(cat1), length(cat2));
    
    for i = 1:length(cat1)
        idx = find(data(:,1) == cat1(i));
        cat2Idx = match(cat2, data(idx,2));
        for j = 1:length(cat2)
            counts(i,j) = length(find(cat2Idx == j));
        end
    end
   
    if PLOTS
        figure, image(counts);

        try
            colormap(jet(maxC));
        catch
            colormap('default');
        end
        %    axis equal;
        title(sprintf('Contingency Table - Non-zero: %2.2f\%', ...
                      100-length(find(counts==0))/(length(cat1)*length(cat2))*100));
        xlabel(text2);
        ylabel(text1);
    end