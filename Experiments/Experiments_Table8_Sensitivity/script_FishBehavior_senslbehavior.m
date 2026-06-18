%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

% parpool(16)
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code
addpath(genpath('.'))
%%
modelPerformance = zeros(15,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
lbhvSetting = [0.5 0.5 0.5 1 1 1 1.01 1.01 1.01 1.1 1.1 1.1 2 2 2];
disp('new 1')
for lbhvSettingIdx = 1:15
    modelPerformance(lbhvSettingIdx,1) = lbhvSetting(lbhvSettingIdx);
    [modelPerformance(lbhvSettingIdx,2), modelPerformance(lbhvSettingIdx,3),modelPerformance(lbhvSettingIdx,4)] = ...
        runBDLDSonZebrafishlbhv(lbhvSetting(lbhvSettingIdx));
end
save('SensitivityExperiment_FishBehavior_senslbehavior.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [medianDynCoefUseOneRun,varExplBhv,varExpl] = runBDLDSonZebrafishlbhv(lambda_behavior)

ifForCIS = 1;

if ifForCIS
    addpath(genpath('.'))
    load("saveFish_En_250901_lbhv1_5000i.mat");
    
else
    addpath(genpath('.'))
    load("saveFish_En_250901_lbhv1_5000i.mat");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inf_opts.lambda_D = 0;

inf_opts.lambda_behavior = lambda_behavior;
inf_opts.max_iters = 500;

inf_opts

%%
% if ifForCIS
%     addpath(genpath('~/my_documents/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code'))
% end
tic

currentState = rng('shuffle');

if inf_opts.behaviordLDS
    [Phi, F, Psi] = bpdndf_dynamics_learning_behavior(dFF, [], [], [], behavior_data, inf_opts);              % Run the dictionary learning algorithm
else
    [Phi, F] = bpdndf_dynamics_learning(dFF, [], [], inf_opts);              % Run the dictionary learning algorithm
end
%%
% if inf_opts.behaviordLDS
%     if inf_opts.AcrossIndividuals
%         for ii = 1:size(Phi,1)
%             [A_cell{ii},B_cell{ii}] = parallel_bilinear_dynamic_inference(dFF(ii), Phi{ii}, F, ...
%                                            @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
%         end
%     else
%         [A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFF, Phi, F, ...
%                                            @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
% 
%     end
if inf_opts.behaviordLDS
    if inf_opts.AcrossIndividuals % EY added 08/27/2023
        for ii = 1:size(data_obj,1)
            [A_cell{ii},B_cell{ii}] = parallel_bilinear_dynamic_inference_behavior(dFF(ii), Phi{ii}, F, Psi{ii}, behaviorData{ii}, ...
                                       @bpdndf_bilinear_handle_behavior, inf_opts); % Infer sparse coefficients
        end
    else
        [A_cell,B_cell] = parallel_bilinear_dynamic_inference_behavior(dFF, Phi, F, Psi, behavior_data, ...
                                       @bpdndf_bilinear_handle_behavior, inf_opts); % Infer sparse coefficients
    end
else
    if inf_opts.AcrossIndividuals
        for ii = 1:size(Phi,1)
            [A_cell{ii},B_cell{ii}] = parallel_bilinear_dynamic_inference(dFF(ii), Phi{ii}, F, ...
                                           @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
        end
    else
        [A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFF, Phi, F, ...
                                           @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
    
    end
end

dLDSruntime = toc;

fprintf('Elapsed time: %.2f seconds\n', dLDSruntime);


%%
varExpl = overallVarExpl(dFF,A_cell,Phi,inf_opts);

%varExplPerNeuron = perNeuronVarExpl(dFF,A_cell,Phi,inf_opts);


medianDynCoefUseOneRun = median(sum(abs(B_cell{1})>1e-3));

reconstructedBehavior = Psi*B_cell{1};
allbhv = behavior_data{1};
rval = corrcoef(allbhv(:), reconstructedBehavior(:));
varExplBhv = (rval(1,2))^2;


end