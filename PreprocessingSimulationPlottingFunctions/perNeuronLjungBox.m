function [lbPerNeuronOriginalS,lbPerNeuronReconstrS,lbPerNeuronMinusReconstrS] = perNeuronLjungBox(dFF,A_cell,Phi,inf_opts,whichCoords,varargin)
if nargin > 5
    lagsRange = varargin{1};
else
    lagsRange = [1 20];
end
fprintf('\nLags %d to %d',lagsRange(1),lagsRange(2));

if inf_opts.AcrossIndividuals
    disp('These dims can be different for every individual - must'+...+
        'be cell')
else
    lbPerNeuronOriginal        = zeros(size(dFF{1},1),size(A_cell,1)); % neurons by trials
    lbPerNeuronReconstr        = zeros(size(dFF{1},1),size(A_cell,1)); % neurons by trials
    lbPerNeuronMinusReconstr   = zeros(size(dFF{1},1),size(A_cell,1)); % neurons by trials

    lbPerNeuronOriginalP       = zeros(size(dFF{1},1),size(A_cell,1)); % neurons by trials
    lbPerNeuronReconstrP       = zeros(size(dFF{1},1),size(A_cell,1)); % neurons by trials
    lbPerNeuronMinusReconstrP  = zeros(size(dFF{1},1),size(A_cell,1)); % neurons by trials

    lbPerNeuronOriginalS       = zeros(size(dFF{1},1),size(A_cell,1)); % neurons by trials
    lbPerNeuronReconstrS       = zeros(size(dFF{1},1),size(A_cell,1)); % neurons by trials
    lbPerNeuronMinusReconstrS  = zeros(size(dFF{1},1),size(A_cell,1)); % neurons by trials

    for ii = 1:size(A_cell,1) %in the event of multiple trials
        disp(ii)
        x_values = A_cell{ii};

        singleDFF = dFF{ii}; % save RAM
        singleReconstr = (Phi * x_values);

        a = double(singleDFF);
        b = double(singleReconstr);
        for jj = 1:size(a,1)
            % figure
            % tiledlayout(2,1)
            % nexttile
            % autocorr(a(jj,:))
            % nexttile
            % parcorr(a(jj,:))
            % 
            % figure
            % tiledlayout(2,1)
            % nexttile
            % autocorr(b(jj,:))
            % nexttile
            % parcorr(b(jj,:))

            [lbPerNeuronOriginal(jj,ii),lbPerNeuronOriginalP(jj,ii)...
                lbPerNeuronOriginalS(jj,ii),~]...
                = lbqtest1d(a(jj,:),'Lags',lagsRange);
            [lbPerNeuronReconstr(jj,ii),lbPerNeuronReconstrP(jj,ii),...
                lbPerNeuronReconstrS(jj,ii),~]...
                = lbqtest1d(b(jj,:),'Lags',lagsRange);
            [lbPerNeuronMinusReconstr(jj,ii),lbPerNeuronMinusReconstrP(jj,ii),...
                lbPerNeuronMinusReconstrS(jj,ii),~]...
                = lbqtest1d(a(jj,:)-b(jj,:),'Lags',lagsRange);

            % aZscored = (a(jj,:)-mean(a(jj,:)))/std(a(jj,:));
            % bZscored = (b(jj,:)-mean(b(jj,:)))/std(b(jj,:));
            % 
            % figure
            % tiledlayout(2,1)
            % nexttile
            % autocorr(bZscored)
            % nexttile
            % parcorr(bZscored)

            % [lbPerNeuronOriginal(jj,ii),lbPerNeuronOriginalP(jj,ii),~,~]...
            %     = lbqtest(aZscored);
            % [lbPerNeuronReconstr(jj,ii),lbPerNeuronReconstrP(jj,ii),~,~]...
            %     = lbqtest(bZscored);
            % [lbPerNeuronVsReconstr(jj,ii),lbPerNeuronVsReconstrP(jj,ii),...
            %     ~,~]...
            %     = lbqtest(aZscored-bZscored);

            close all
        end

        figure()
        subplot(3,3,1)
        histogram(lbPerNeuronOriginal(:),7)
        title('h vals, Original Data')
        subplot(3,3,2)
        histogram(lbPerNeuronOriginalP(:),15)
        title('p vals, Original Data')
        subplot(3,3,3)
        histogram(lbPerNeuronOriginalS(:),15)
        title('Test stat, Original Data')
        subplot(3,3,4)
        histogram(lbPerNeuronReconstr(:),7)
        title('h vals, Reconstructed Data')
        subplot(3,3,5)
        histogram(lbPerNeuronReconstrP(:),15)
        title('p vals, Reconstructed Data')
        subplot(3,3,6)
        histogram(lbPerNeuronReconstrS(:),15)
        title('Test stat, Reconstructed Data')
        subplot(3,3,7)
        histogram(lbPerNeuronMinusReconstr(:),7)
        title('h vals, Orig - Reconstr Data')
        subplot(3,3,8)
        histogram(lbPerNeuronMinusReconstrP(:),15)
        title('p vals, Orig - Reconstr Data')
        subplot(3,3,9)
        histogram(lbPerNeuronMinusReconstrS(:),15)
        title('Test stat, Orig - Reconstr Data')
        set(gcf,"Color","white")

        figure()
        hold on
        colormap(parula)
        scatter3(whichCoords(:,3),whichCoords(:,2),whichCoords(:,1),50, lbPerNeuronReconstrS(:)-lbPerNeuronOriginalS(:))
        colorbar
        title('Ljung-Box delta(test statistics) per neuron')
        set(gcf,"Color","white")
        hold off
    end
end

end