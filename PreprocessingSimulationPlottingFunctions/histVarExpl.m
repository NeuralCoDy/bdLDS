function histVarExpl(whichData,inf_opts,trialIdentifier,varargin)
if nargin > 3
    dataIdentifier = varargin{1};
else
    dataIdentifier = "zebgen";
end

if inf_opts.AcrossIndividuals
    for i = 1:size(whichData,1) % individual 
        fig1 = figure();
        cla
        hold on
        histogram(whichData{i}(:))
        title('Var expl per neuron')
        set(gcf,"Color","white")
        hold off
        filename = sprintf("%s_%s_histvarexpl_%d.fig", dataIdentifier, trialIdentifier,i);
        saveas(fig1,filename);
        %close all
    end

else

    for i = 1:size(whichData,2) % trial
        fig1 = figure();
        cla
        hold on
        histogram(whichData(:,i))
        title('Var expl per neuron')
        set(gcf,"Color","white")
        hold off
        filename = sprintf("%s_%s_histvarexpl_%d.fig", dataIdentifier, trialIdentifier,i);
        saveas(fig1,filename);
        %close all
    end
end

end