function varExpl = varExplLatentsOnly(inferredLatents, groundTruthLatents)

if size(inferredLatents{1},1) ~= size(groundTruthLatents{1},1)
    varExpl = cell(size(inferredLatents,1),1); %,size(inferredLatents{1},1),size(groundTruthLatents{1},1)
else
    varExpl = zeros(1,size(inferredLatents,1));
end


for ii = 1:size(inferredLatents,1) %in the event of multiple trials
    if size(inferredLatents{ii},1) ~= size(groundTruthLatents{ii},1)
        for jj = 1:size(inferredLatents{ii},1) 
            for kk = 1:size(groundTruthLatents{ii},1)
                rval = corrcoef(inferredLatents{ii}(jj,:), groundTruthLatents{ii}(kk,:));
                varExpl{ii}(jj,kk) = (rval(1,2))^2;
            end
        end
    else
        rval = corrcoef(inferredLatents{ii}, groundTruthLatents{ii});
        varExpl(ii) = (rval(1,2))^2;
    end
    % pause();
end
disp(varExpl);

end