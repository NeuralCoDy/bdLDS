%%
load("SensitivityExperiment_Sim_lambdaDynCoef.mat")

ifZebrafish = 0; % change if showing NoObs sim data, since x reconstruction is irrelevant
%%
figure();
x0=10;
y0=10;
width=650;
height=250;
set(gcf,'units','points','position',[x0,y0,width,height])
tiledlayout("horizontal")
modelPerformance_var01 = reshape(modelPerformance(:,1),3,5); % param sweep settings
modelPerformance_var02 = reshape(modelPerformance(:,2),3,5); % median dyn coef use
modelPerformance_var03 = reshape(modelPerformance(:,3),3,5); % behavior R^2

paramSweepVarName = "lambda dyn coef";

if ifZebrafish
    modelPerformance_var04 = reshape(modelPerformance(:,4),3,5); % neural R^2
    nexttile()
    %errorbar(modelPerformance_var01(1,:),mean(modelPerformance_var02),std(modelPerformance_var02)/sqrt(3),'DisplayName','b-dLDS','LineStyle', 'none', 'Marker', 'o')
    errorbar(modelPerformance_var01(1,:),mean(modelPerformance_var02),std(modelPerformance_var02)/sqrt(3),"-o","MarkerSize",5,"MarkerFaceColor",[1 1 1],'DisplayName','b-dLDS','LineStyle', '--')
    ylabel('median dyn coef use')
    xlabel(paramSweepVarName)
    box off
    nexttile()
    errorbar(modelPerformance_var01(1,:),mean(modelPerformance_var03),std(modelPerformance_var03)/sqrt(3),"-o","MarkerSize",5,"MarkerFaceColor",[1 1 1],'DisplayName','b-dLDS','LineStyle', '--')
    ylabel('behavior R^2')
    xlabel(paramSweepVarName)
    box off
    nexttile()
    errorbar(modelPerformance_var01(1,:),mean(modelPerformance_var04),std(modelPerformance_var04)/sqrt(3),"-o","MarkerSize",5,"MarkerFaceColor",[1 1 1],'DisplayName','b-dLDS','LineStyle', '--')
    ylabel('neural R^2')
    xlabel(paramSweepVarName)
    box off
else
    disp('If it looks like points have disappeared, check for NaNs; omitmissing not used')
    nexttile()
    errorbar(modelPerformance_var01(1,:),mean(modelPerformance_var02),std(modelPerformance_var02)/sqrt(3),"-o","MarkerSize",5,"MarkerFaceColor",[1 1 1],'DisplayName','b-dLDS','LineStyle', '--')
    ylabel('median dyn coef use')
    xlabel(paramSweepVarName)
    box off
    nexttile()
    errorbar(modelPerformance_var01(1,:),mean(modelPerformance_var03),std(modelPerformance_var03)/sqrt(3),"-o","MarkerSize",5,"MarkerFaceColor",[1 1 1],'DisplayName','b-dLDS','LineStyle', '--')
    ylabel('behavior R^2')
    xlabel(paramSweepVarName)
    box off

end