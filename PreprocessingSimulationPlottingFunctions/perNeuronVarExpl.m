function varExplPerNeuron = perNeuronVarExpl(dFF,A_cell,Phi,inf_opts)
if inf_opts.AcrossIndividuals
    disp('These dims can be different for every individual - must'+...+
        'be cell')

    % for ii = 1:size(Phi,1)
	%     varExplPerNeuron{ii}     = zeros(size(dFF{ii},1),1);
    %     x_values = cell2mat(A_cell{ii});
    %     dataReconstructed = (Phi{ii} * x_values); 
	%     disp(size(dFF))
	%     disp(size(dataReconstructed))
    % 
    %     singleDFF = single(dFF{ii}(1:10,:)); % save RAM
    %     singleReconstr = single(dataReconstructed(1:10,:));
    % 
    %     a = singleDFF;
    %     b = singleReconstr;
    % 
    %     disp('check if dims are correct')
    %     varExplPerNeuron{ii}(:) = diag(((1/(size(a,2)-1)*...
    %         ( ((a-mean(a,2))./std(a,0,2))*((b-mean(b,2))./std(b,0,2)).')).^2));
    % end
else
    varExplPerNeuron     = zeros(size(dFF{1},1),size(A_cell,1)); % neurons by trials

    for ii = 1:size(A_cell,1) %in the event of multiple trials
        disp(ii)
        x_values = A_cell{ii};

        singleDFF = dFF{ii}; % save RAM % N by T
        singleReconstr = (Phi * x_values); % N by T

        a = singleDFF;
        b = singleReconstr;
        varExplPerNeuron(:,ii) = diag(((1/(size(a,2)-1)*...
            ( ((a-mean(a,2))./std(a,0,2))*((b-mean(b,2))./std(b,0,2)).')).^2));

    end
end

end
