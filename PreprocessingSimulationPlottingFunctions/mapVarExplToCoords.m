function mapVarExplToCoords(whichData,whichCoords,inf_opts,trialIdentifier)

if inf_opts.AcrossIndividuals
    for i = 1:size(whichData,1) % individual 
        fig1 = figure();
        cla
        hold on
        colormap(parula)
        scatter3(whichCoords(:,3),whichCoords(:,2),whichCoords(:,1),50, whichData{i}(:))
        colorbar
        caxis([0 1])
        title('Var expl per neuron')
        set(gcf,"Color","white")
        hold off
        filename = sprintf("zebgen_%s_mapvarexpl_%d.fig", trialIdentifier,i);
        saveas(fig1,filename);
        % close all
    end

else

    for i = 1:size(whichData,2) % trial
        fig1 = figure();
        cla
        hold on
        colormap(parula)
        scatter3(whichCoords(:,3),whichCoords(:,2),whichCoords(:,1),50, whichData(:,i))
        colorbar
        %caxis([-2 2])
        caxis([0 1])
        title('Var expl per neuron')
        set(gcf,"Color","white")
        hold off
        filename = sprintf("zebgen_%s_mapvarexpl_%d.fig", trialIdentifier,i);
        saveas(fig1,filename);
        % close all
    end
end

end
