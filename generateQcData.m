function qcSpan = generateQcData(QC, lat, lon, loadTime, truckId, trucks)
% GENERATEQCDATA - Function to generate daily QC data
%   


% Get QC data
QCs = unique(QC);
qcSpan(length(QCs)).latSpan = [];
qcSpan(length(QCs)).lonSpan = [];
qcSpan(length(QCs)).latP = [];
qcSpan(length(QCs)).lonP = [];
qcSpan(length(QCs)).day = [];
qcSpan(length(QCs)).numTickets = [];

disp('Generating QC Data...');

for i = 1:length(QCs)
    loadIdx = find(QC == QCs(i) & lat~=0 & lon~=0);
    
    % Find daily data for each QC
    k = 1;
    for j=floor(min(loadTime)):floor(max(loadTime))
        idx = loadIdx(find(floor(loadTime(loadIdx))==j));
        
        if ~isempty(idx)
            qcSpan(i).day(k) = j;
            qcSpan(i).numTickets(k) = length(idx);
            qcSpan(i).latSpan(k) = max(lat(idx)) - ...
                min(lat(idx)) + 0.0002;
            qcSpan(i).lonSpan(k) = max(lon(idx)) - ...
                min(lon(idx)) + 0.0002;
            qcSpan(i).latP(k) = min(lat(idx)) - 0.0001;
            qcSpan(i).lonP(k) = min(lon(idx)) - 0.0001;
            qcSpan(i).trucks(:,k) = contingencyTable([truckId(idx), ...
                                QC(idx)], trucks, QCs(i), 'a', ...
                                                     'b', 10, 0);
            
            % Increment day
            k = k + 1;
        end
    end
end