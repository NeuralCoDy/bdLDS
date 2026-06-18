function [filteredData, varargout] = filterByStd(dataToFilter,...
    thresholdNumberStd,wholeDataOrByNeuron,zeroOrMeanOrMode,varargin)
disp('Note: filters and caps')
originalData = dataToFilter;
filteredIndices = (1:size(originalData,1)).';

disp('Data loaded')
fprintf('\n Size of dataToFilter: %d %d',size(dataToFilter,1),size(dataToFilter,2))

if nargout > 1
    relatedDataToFilter = varargin{1};
    fprintf('\n Size of relatedDataToFilter: %d %d',size(relatedDataToFilter,1),size(relatedDataToFilter,2))
end
if nargout > 2
    relatedDataToFilter2 = varargin{2};
    fprintf('\n Size of relatedDataToFilter: %d %d',size(relatedDataToFilter,1),size(relatedDataToFilter,2))
end

% deal with NaNs
[rowsToExclude, ~] = find(isnan(dataToFilter));
uniqueRowsToExclude = unique(rowsToExclude);
% disp(uniqueRowsToExclude)
dataToFilter(uniqueRowsToExclude,:) = [];
% dFF(isnan(dFF)) = 0;

if nargout > 1
    relatedDataToFilter(uniqueRowsToExclude,:) = [];
end
if nargout > 2
    relatedDataToFilter2(uniqueRowsToExclude,:) = [];
end

filteredIndices(uniqueRowsToExclude,:) = [];

% thresholdReasonable = 0.1; %reasonable dFF values for zebrafish
% isAllAbove = all( dataToFilter > thresholdReasonable, 2 ) ;
% dataToFilter(isAllAbove,:) = [];
% 
% if nargout > 1
%     relatedDataToFilter(isAllAbove,:) = [];
% end
% if nargout > 2
%     relatedDataToFilter2(isAllAbove,:) = [];
% end


% Goal: filter out rows where all of the values are near "zero"/baseline
if zeroOrMeanOrMode == 0 % zero
    avgThr = 0;
    if wholeDataOrByNeuron
        stdThr = std(dataToFilter(:));
    else
        stdThr = std(dataToFilter,0,2);
    end

elseif zeroOrMeanOrMode == 2 % mode
    if wholeDataOrByNeuron
        modeEst = kernelModeEstimate(dataToFilter(:));
        stdThr = std(dataToFilter(:));
        avgThr = modeEst;
    else
        modeEst = kernelModeEstimate(dataToFilter.');
        stdThr = std(dataToFilter.').';
        avgThr = modeEst.';
    end
    
else % mean
    if wholeDataOrByNeuron
        stdThr = std(dataToFilter(:));
        avgThr = mean(dataToFilter(:));
    else
        stdThr = std(dataToFilter,0,2);
        avgThr = mean(dataToFilter,2);
    end

end

thresholdNearAvgUpper = avgThr + 0.01*stdThr;
thresholdNearAvgLower = avgThr - 0.01*stdThr;

%thresholdToKeepDownUpper = 0;
%thresholdToKeepDownLower = 0;
%thresholdToKeepDown = meanThr - thresholdNumberStd*stdThr;

% Keep a row if at least one value is stronger than near-0 range
[rowsToKeep, ~] = find((dataToFilter > thresholdNearAvgUpper)| (dataToFilter < thresholdNearAvgLower));
uniqueRowsToKeep = unique(rowsToKeep);
% disp(uniqueRowsToKeep)
dataToFilter = dataToFilter(uniqueRowsToKeep,:);
% disp(size(dataToFilter))
if nargout > 1
    relatedDataToFilter = relatedDataToFilter(uniqueRowsToKeep,:);
end
if nargout > 2 % Important: fixed from 3 - EY 12/30/25. Does not affect En's plots because those only had dFF + coordinates, but does affect Emmanuel's plots (data, genes, coords). As of edits, those were not finalized yet anyway.
    relatedDataToFilter2 = relatedDataToFilter2(uniqueRowsToKeep,:);
end

filteredIndices = filteredIndices(uniqueRowsToKeep,:);

fprintf('\n Size of dataToFilter: %d %d',size(dataToFilter,1),size(dataToFilter,2))
if nargout > 1
    fprintf('\n Size of relatedDataToFilter: %d %d',size(relatedDataToFilter,1),size(relatedDataToFilter,2))
end
if nargout > 2 
    fprintf('\n Size of relatedDataToFilter2: %d %d',size(relatedDataToFilter2,1),size(relatedDataToFilter2,2))
end

% Set upper and lower thresholds for extreme values that should be capped
if size(avgThr,1) > 1
    avgThr = avgThr(uniqueRowsToKeep,1);
end

if size(stdThr,1) > 1
    stdThr = stdThr(uniqueRowsToKeep,1);
end

thresholdExtremeUpper = ones(size(uniqueRowsToKeep,1),1);%avgThr + 5*stdThr;
thresholdExtremeLower = -1*ones(size(uniqueRowsToKeep,1),1);%avgThr - 5*stdThr;

% % all values not too extreme (large positive or negative)
% % The trouble is that this removes all but 56 records - nearly all of the 
% % traces contain at least one extreme value. Moreover, the extreme values
% % are not measurement problems - they are technically right and possible.
% [rowsToExclude, ~] = find(dataToFilter>thresholdExtremeUpper);
% uniqueRowsToExclude = unique(rowsToExclude);
% dataToFilter(uniqueRowsToExclude,:) = [];
% if nargin > 4
%     relatedDataToFilter(uniqueRowsToExclude,:) = [];
% end
% if nargin > 5
%     relatedDataToFilter2(uniqueRowsToExclude,:) = [];
% end
% 
% [rowsToExclude, ~] = find(dataToFilter<thresholdExtremeLower);
% uniqueRowsToExclude = unique(rowsToExclude);
% dataToFilter(uniqueRowsToExclude,:) = [];
% if nargin > 4
%     relatedDataToFilter(uniqueRowsToExclude,:) = [];
% end
% if nargin > 5
%     relatedDataToFilter2(uniqueRowsToExclude,:) = [];
% end

% cap individual values at -0.1 and 0.1
[rowsToCap,~] = find(dataToFilter>thresholdExtremeUpper);
idxToCap      = find(dataToFilter>thresholdExtremeUpper);
dataToFilter(idxToCap) = thresholdExtremeUpper(rowsToCap);

[rowsToCap,~] = find(dataToFilter<thresholdExtremeLower);
idxToCap      = find(dataToFilter<thresholdExtremeLower);
dataToFilter(idxToCap) = thresholdExtremeLower(rowsToCap);

disp('After capping')
fprintf('\n Size of dataToFilter: %d %d',size(dataToFilter,1),size(dataToFilter,2))
if nargout > 1
    fprintf('\n Size of relatedDataToFilter: %d %d',size(relatedDataToFilter,1),size(relatedDataToFilter,2))
end
if nargout > 2
    fprintf('\n Size of relatedDataToFilter2: %d %d',size(relatedDataToFilter2,1),size(relatedDataToFilter2,2))
end

% filter out rows with std = 0
disp('new: filtering by std Nov 6 2025')
[rowsConstant,~] = find(std(dataToFilter,0,2)==0);
disp(rowsConstant)
disp('std0')
disp(size(rowsConstant))
dataToFilter(rowsConstant,:) = [];

stdDFF2 = robustSTD(dataToFilter.');
[rowsConstant2,~] = find(stdDFF2==0);
disp(rowsConstant2)
disp('robust std')
disp(size(rowsConstant2))
dataToFilter(rowsConstant2,:) = [];

% IMPORTANT: I missed this before! Added 12/30/25
if nargout > 1
    relatedDataToFilter(rowsConstant,:)  = [];
    relatedDataToFilter(rowsConstant2,:)  = [];
end
if nargout > 2
    relatedDataToFilter2(rowsConstant,:) = [];
    relatedDataToFilter2(rowsConstant2,:) = [];
end
filteredIndices(rowsConstant,:)          = [];
filteredIndices(rowsConstant2,:)          = [];


% check 
max(dataToFilter,[],"all")
min(dataToFilter,[],"all")

% 
% 
% if nargin > 4
%     varargout{1} = relatedDataToFilter(uniqueRowsToKeep,:);
% end
% if nargin > 5
%     varargout{2} = relatedDataToFilter2(uniqueRowsToKeep,:);
% end

disp('After filtering by std')
fprintf('\n Size of dataToFilter: %d %d',size(dataToFilter,1),size(dataToFilter,2))
if nargout > 1
    fprintf('\n Size of relatedDataToFilter: %d %d',size(relatedDataToFilter,1),size(relatedDataToFilter,2))
end
if nargout > 2
    fprintf('\n Size of relatedDataToFilter2: %d %d',size(relatedDataToFilter2,1),size(relatedDataToFilter2,2))
end


filteredData = dataToFilter;
disp(size(filteredData))
if nargout > 1
    varargout{1} = relatedDataToFilter;
end
if nargout > 2
    varargout{2} = relatedDataToFilter2;
end
if nargout > 3
    varargout{3} = filteredIndices;
end

end
