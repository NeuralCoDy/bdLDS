%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

% parpool(16)
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code
addpath(genpath('.'))
%%
modelPerformance = zeros(15,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
lambdaDynCoefSetting = [0.25 0.25 0.25 0.251 0.251 0.251 0.26 0.26 0.26 0.35 0.35 0.35 1.35 1.35 1.35];
for lambdaDynCoefSettingIdx = 1:15
    modelPerformance(lambdaDynCoefSettingIdx,1) = lambdaDynCoefSetting(lambdaDynCoefSettingIdx);
    [modelPerformance(lambdaDynCoefSettingIdx,2), modelPerformance(lambdaDynCoefSettingIdx,3), modelPerformance(lambdaDynCoefSettingIdx,4)] = ...
        runBDLDSonSim(lambdaDynCoefSetting(lambdaDynCoefSettingIdx));
end

save('SensitivityExperiment_Sim_lambdaDynCoef.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [medianDynCoefUseOneRun,varExplBhv,varExpl] = runBDLDSonSim(lambda_b)
nBhv = 1;
howManyColsPhi = 1;
noiseLevel = 0;
inf_opts.lambda_b = lambda_b;

%%
rng('shuffle')
currentState = rng;
%%
sizeD            = 8;
sizeT            = 3000; % T/2
simulatedDmatrix = eye(sizeD,sizeD);
contOpt          = false;
nF               = [3,3];

[latentStatesX, Fgt, groundTruthStates] = generateContMultiSubNetwork_Range(...
               sizeT,[sizeD/2,sizeD/2],nF,50,contOpt,false,false,true);
%%
% nBhv = 1;
A = randn(nBhv,sum(nF,'all')); 
if howManyColsPhi < 6
    A(:,howManyColsPhi+1:end) = 0; 
end
simulatedPsi = A * 1./vecnorm(A);
simulatedPsi(isnan(simulatedPsi))=0;

% sum(simulatedPsi.^2,1)
% simulatedPsi.'*simulatedPsi
% corr(simulatedPsi)
% disp('check for correlation between columns - if too high, might impede reconstruction of behavior')
justOneGTState = groundTruthStates;
behaviorData = cellfun(@(x) simulatedPsi*x,justOneGTState,'UniformOutput',false);
%%
% noiseLevel = input('noise level (default:0):');
% if isempty(noiseLevel)
%     noiseLevel = 0;
% end
dFF = statesToObservations(simulatedDmatrix, latentStatesX, noiseLevel);

%%
% Set parameters
inf_opts.nF              = 15; % # dynamics operators
inf_opts.N               = size(latentStatesX{1},1); % # latent states
inf_opts.M               = size(dFF{1},1) ; % original dimension (# channels)
% disp('If reset dFF below, make sure to change inf_opts.M');
inf_opts.lambda_val      = 0.0001;         % ASC added 7/25 % VARY - tradeoff: approx. SNR - see ratios of traces % 0.0001 (don't shrink this further - gets to be essentially 0 for solver_L1RLS)
inf_opts.lambda_history  = 0.0001;%0.01;         % ASC added 7/25
% inf_opts.lambda_b        = 0.25; %0.5;%0.39;           % ASC added 7/25 %0.025 % 0.5 %0.6
inf_opts.lambda_historyb = 0.45;%1.16;           % ASC added 7/25 %0.7 %0.45 %0.35 % 1.99 way too strong - all coefficients end up same
inf_opts.tol             = 1e-8;           % 1e-3
inf_opts.max_iter2       = 20;             % 500 %20
inf_opts.max_iters       = 500; %input('max iters:');
fprintf('\nmax_iters: %d\n',inf_opts.max_iters)
inf_opts.F_update        = true; % default = true;
inf_opts.D_update        = false; % NoObs case - default = true;
inf_opts.N_ex            = 100; %input('N_ex (e.g., 50):');%50;   
inf_opts.T_s             = 200; %input('T_s (e.g., 200):');%200; 
if inf_opts.T_s > size(dFF{1},2)
    warning('Sample duration (time points T_s) too large for data')
    fprintf('\nDefaulting to sample size == size of data %d\n',size(dFF{1},2))
    inf_opts.T_s = size(dFF{1},2);
end
inf_opts.step_d          = 1;  % ASC added 7/25
inf_opts.step_f          = 10;  % ASC added 7/25 % 30
inf_opts.step_decay      = 0.99995; %input('step_decay (e.g., 0.9999):'); %0.998 %0.9999 helped get closer - 4/6 operators in corners 
inf_opts.plot_option     = 10;  % ASC added 7/25
inf_opts.lambda_f_decay  = 0.996; %0.996
inf_opts.lambda_f        = 0.1; %input('lambda_f (e.g., 0.2):'); %0.2 almost there but two operators nearly identical %0.4 too high - only 1-2 kinds of operators emerge
inf_opts.solver_type     = 'fista'; %input('Solver (fista or tfocs)?','s');%'tfocs'; % fista
inf_opts.special         = 'noobs';
inf_opts.deltaDynamics   = false; %default: false %use x_t-x_{t-1}
inf_opts.AcrossIndividuals     = 0; %input('Across individuals? 1 for yes, 0 for no:');

inf_opts.behaviordLDS    = 1; %input('behavior_dLDS?:');
if inf_opts.behaviordLDS == 1
    inf_opts.step_psi        = 10; %input('step_psi:');%1; 
    inf_opts.lambda_behavior = 0.1; %input('lambda_behavior:'); %1.6; %0.5521; % 0.5536 for opt for psi
    inf_opts.verysparsebhv   = 1;
    inf_opts.lambda2         = 5; %input('lambda2:');
    inf_opts.psinorm         = 'frob'; %input('Norm on Psi - norm or frob?:','s');
else
    inf_opts.step_psi        = 0; % 30 for opt for psi
    inf_opts.lambda_behavior = 0; %0.5521; % 0.5536 for opt for psi
    inf_opts.verysparsebhv   = 0;
end

inf_opts

%% because CLDS needs an 80/20 train/test split
whichTrials   = randsample(50,40);
dFF2          = dFF(whichTrials);
behaviorData2 = behaviorData(whichTrials);
%%
if inf_opts.behaviordLDS
    [Phi, F, Psi] = bpdndf_dynamics_learning_behavior(dFF2, [], [], [], behaviorData2, inf_opts);              % Run the dictionary learning algorithm
else
    [Phi, F] = bpdndf_dynamics_learning(dFF2, [], [], inf_opts);              % Run the dictionary learning algorithm
end
%%
if inf_opts.behaviordLDS
    if inf_opts.AcrossIndividuals % EY added 08/27/2023
        for ii = 1:size(data_obj,1)
            [A_cell{ii},B_cell{ii}] = parallel_bilinear_dynamic_inference_behavior(dFF(ii), Phi{ii}, F, Psi{ii}, behaviorData{ii}, ...
                                       @bpdndf_bilinear_handle_behavior, inf_opts); % Infer sparse coefficients
        end
    else
        [A_cell,B_cell] = parallel_bilinear_dynamic_inference_behavior(dFF, Phi, F, Psi, behaviorData, ...
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

medianDynCoefUseOneRun = median(sum(abs(B_cell{1})>1e-3));

reconstructedBehavior = Psi*B_cell{1};
allbhv = behaviorData{1};
rval = corrcoef(allbhv(:), reconstructedBehavior(:));
varExplBhv = (rval(1,2))^2;

varExpl = overallVarExpl(dFF,A_cell,Phi,inf_opts);

end