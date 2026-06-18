function varExpl = overallVarExpl(dFF,A_cell,Phi,inf_opts)

if inf_opts.AcrossIndividuals
    for ii = 1:size(Phi,1)
        x_values = cell2mat(A_cell{ii});
        dataReconstructed = (Phi{ii} * x_values).';    
        
        flatDFF = reshape(dFF{ii}, numel(dFF{ii}),[]);
        flatReconstructed = reshape(dataReconstructed,numel(dataReconstructed),[]);
        rval = corrcoef(dFF{ii}, dataReconstructed.');
        varExpl = (rval(1,2))^2;
        disp(varExpl);

    end
else
    for ii = 1:size(A_cell,1) %in the event of multiple trials
        disp(ii)
        x_values = A_cell{ii};
        %dataReconstructed = (Phi * x_values).';    
        
        maxvalDFF = max(dFF{ii},[],'all');
        stdzdOn = false;
        if stdzdOn
            dataReconstructedRescaled = maxvalDFF .*(Phi * x_values);
            dFFRescaled = maxvalDFF .* dFF{ii};
        else
            dFFRescaled = dFF{ii};
            dataReconstructedRescaled = (Phi * x_values);
        end
       
        
        flatDFF = reshape(dFFRescaled, numel(dFFRescaled),[]);
        flatReconstructed = reshape(dataReconstructedRescaled,numel(dataReconstructedRescaled),[]);
        rval = corrcoef(flatDFF, flatReconstructed);
        varExpl = (rval(1,2))^2;
        disp(varExpl);
        % pause();
    end
end


end