%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% List of tasks:
% Rate estimation?
% Basic dimensionality reduction?
% Dynamics estimation?

% parpool(2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code
addpath(genpath('.'))
addpath(genpath('./External_Packages/distinguishable_colors/'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%0000000l;%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up pointers to the fluorescence data and meta behavioral data

% dat.dir   = 'C:/Users/Helen/Documents/Documents/Documents (3)/GradSchool/JohnsHopkinsStudent/CharlesLab/CElegans/Data';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data (as a cell array)

close all;
%%
clear all;
%%
todaysDate = input('Date (e.g., 251126):','s');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% addpath(genpath('~/my_documents/dLDS-Discrete-Matlab-Model_1FishManyTrials/code'))
%%
% disp('Change path to modeEstimates for your setup')
% addpath(genpath('~/my_documents/CoDybase/CoDybase-MATLAB/CoDybase-MATLAB/stats/'))
ifForCIS = input('CIS?(1 or 0):');
if ifForCIS
    addpath(genpath('~/my_documents/CoDybase/CoDybase-MATLAB/CoDybase-MATLAB/'))
else
    addpath(genpath('C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/CoDybase/CoDybase-MATLAB/CoDybase-MATLAB/'))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up pointers to the fluorescence data and meta behavioral data

%dat.dir = '../../../adamsc/data/zebrafish/';
if ifForCIS
    addpath(genpath('.'))
    load("saveFish_En_250901_lbhv1_5000i.mat");
    % % Add all the required packages that are custom code
    % addpath(genpath('../../npy-matlab/npy-matlab/'))                % Use steinmetz code for loady npy into matlab
    % 
    % % addpath(genpath('../../../../home/adamsc/GITrepos/zebrafishDynamics/code/'))            % This is the main codebase for the project
    % % addpath(genpath('../../../../home/adamsc/GITrepos/dynamics_learning/'))                 % This is the package for learning dynamics 
    % 
    % addpath(genpath('../../../adamsc/data/zebrafish/'))
    % 
    % dat.dir   = './FishBehaviorData/';
    % dat.Ffile = 'for_JHU_cells_dff.npy';
    % dat.Bfile = 'for_JHU_gad1b-6hz-backwardpulse_integrate_behavior_ds.npy';
    % dat.Vfile = 'for_JHU_-gad1b-6hz-backwardpulse_integrate_visual_velocity_ds.npy';
else
    addpath(genpath('.'))
    load("saveFish_En_250901_lbhv1_5000i.mat");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inf_opts.lambda_D = 0;

inf_opts.verysparsebhv = 1;
inf_opts.psinorm       = 'frob';
inf_opts.lambda2       = 5;

inf_opts.lambda_history  = 0.1154; % BayesOpt
inf_opts.lambda_behavior = 1.9953; % BayesOpt

inf_opts.nF = 50;
%%
% if ifForCIS
%     addpath(genpath('~/my_documents/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code'))
% end
tic

% https://www.mathworks.com/matlabcentral/fileexchange/51945-memorylinux
[~,meminfo] = system('cat /proc/meminfo'); 
tokens = regexpi(meminfo,'^MemTotal:\s*(\d+)\s', 'tokens'); 
totalmem = str2double(tokens{1}{1});  
 % get available memory                                                  
tokens = regexpi(meminfo,'^*MemFree:\s*(\d+)\s','tokens');    
freemem = str2double(tokens{1}{1});                 
startmem = totalmem-freemem;    

currentState = rng('shuffle');

if inf_opts.behaviordLDS
    [Phi, F, Psi] = bpdndf_dynamics_learning_behavior_x(dFF, [], [], [], behavior_data, inf_opts);              % Run the dictionary learning algorithm
else
    [Phi, F] = bpdndf_dynamics_learning(dFF, [], [], inf_opts);              % Run the dictionary learning algorithm
end

if inf_opts.behaviordLDS
    if inf_opts.AcrossIndividuals % EY added 08/27/2023
        for ii = 1:size(data_obj,1)
            [A_cell{ii},B_cell{ii}] = parallel_bilinear_dynamic_inference_behavior(dFF(ii), Phi{ii}, F, Psi{ii}, behaviorData{ii}, ...
                                       @bpdndf_bilinear_handle_behavior_x, inf_opts); % Infer sparse coefficients
        end
    else
        [A_cell,B_cell] = parallel_bilinear_dynamic_inference_behavior(dFF, Phi, F, Psi, behavior_data, ...
                                       @bpdndf_bilinear_handle_behavior_x, inf_opts); % Infer sparse coefficients
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

% https://www.mathworks.com/matlabcentral/fileexchange/51945-memorylinux
[~,meminfo] = system('cat /proc/meminfo'); 
tokens = regexpi(meminfo,'^MemTotal:\s*(\d+)\s', 'tokens'); 
totalmem = str2double(tokens{1}{1});  
 % get available memory                                                  
tokens = regexpi(meminfo,'^*MemFree:\s*(\d+)\s','tokens');    
freemem = str2double(tokens{1}{1});                 
endmem = totalmem-freemem; 

memused = endmem-startmem;

varExpl = overallVarExpl(dFF,A_cell,Phi,inf_opts);
%%
varExplPerNeuron = perNeuronVarExpl(dFF,A_cell,Phi,inf_opts);
%%
save(sprintf('saveFish_En_%s_bdLDSx_s07_psireg_nF50.mat',todaysDate),'-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%