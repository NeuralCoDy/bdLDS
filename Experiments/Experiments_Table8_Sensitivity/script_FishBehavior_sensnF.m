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
nFSetting = [5 5 5 10 10 10 20 20 20 25 25 25 50 50 50];
disp('new')
for nFSettingIdx = 1:15
    modelPerformance(nFSettingIdx,1) = nFSetting(nFSettingIdx);
    [modelPerformance(nFSettingIdx,2), modelPerformance(nFSettingIdx,3),modelPerformance(nFSettingIdx,4)] = ...
        runBDLDSonZebrafishnF(nFSetting(nFSettingIdx));
end
save('SensitivityExperiment_FishBehavior_sensnF.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [medianDynCoefUseOneRun,varExplBhv,varExpl] = runBDLDSonZebrafishnF(nF)

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

inf_opts.nF = nF;
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