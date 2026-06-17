function boutType = calculateZebrafishBoutType(startPoint,endPoint,boutData) %92, 5041
nPoints = endPoint-startPoint+1;
boutDuration = 150; %specific to the zebrafish setup
boutType = strings(1,nPoints); %5041-91

nBouts = round(nPoints/boutDuration);

velTrimmed = boutData(startPoint:endPoint); % https://www.mathworks.com/matlabcentral/answers/326675-cumulative-sum-at-over-an-specified-interval
B=reshape(velTrimmed,boutDuration,[]);
D=string(reshape(velTrimmed,boutDuration,[]));

for idxBout = 1:nBouts
    boutTypeThisBout = "";
    appliedFirstDisplacement = sum(B(1:3,idxBout));   %specific to the zebrafish setup
    appliedSecondDisplacement = sum(B(40:54,idxBout));  
    if appliedFirstDisplacement > 0 & appliedSecondDisplacement > 0
        boutTypeThisBout = "+Both";
    elseif appliedFirstDisplacement > 0 & appliedSecondDisplacement < 0
        boutTypeThisBout = "+First-Second (Total 0)";
    elseif appliedFirstDisplacement < 0 & appliedSecondDisplacement < 0
        boutTypeThisBout = "-Both";
    elseif appliedFirstDisplacement > 0 & (appliedSecondDisplacement-0) < 0.0001
        boutTypeThisBout = "+First";
    elseif appliedFirstDisplacement < 0 & (appliedSecondDisplacement-0) < 0.0001
        boutTypeThisBout = "-First";
     elseif appliedSecondDisplacement > 0 & (appliedFirstDisplacement-0) < 0.0001
        boutTypeThisBout = "+Second";
    elseif appliedSecondDisplacement < 0 & (appliedFirstDisplacement-0) < 0.0001
        boutTypeThisBout = "-Second";
    else
        boutTypeThisBout = "No applied displacement";
    end
        

    D(:,idxBout) = repelem(boutTypeThisBout,boutDuration).'; %FIXME
end

boutType=reshape(D,[],1);


end