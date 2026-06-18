function varToFlip = checkAndFlip(varToFlip,timepoints)
        if size(varToFlip,2) == size(timepoints,1)
            varToFlip = varToFlip.'; 
        end
end

