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

inf_opts.lambda_D  = 0;
inf_opts.max_iters = 100;
inf_opts.nF        = 50; % to match inf_opts.N, for comparison of coefficient R^2 vs. behavior
inf_opts.N         = 50;
inf_opts.psinorm   = "frob";
inf_opts.verysparsebhv   = 1;

%% bayesopt

rng default
lambda_val      = optimizableVariable('lambda_val',[0.01 1],'Type','real');
lambda_history  = optimizableVariable('lambda_history',[0.01 1],'Type','real'); % [0.001 1] --> 0.5536
lambda_b        = optimizableVariable('lambda_b',[0.01 1],'Type','real');
lambda_behavior = optimizableVariable('lambda_behavior',[0.01 4],'Type','real'); % [0.1 10], real --> 9.9994
step_psi        = optimizableVariable('step_psi',[1 20],'Type','integer');
lambda2         = optimizableVariable('lambda2',[0.1 20],'Type','real');
fun = @(x)runBehaviorDLDSBO_opt6params_zebrafish(dFF,behavior_data,...
    x.lambda_val,x.lambda_history,x.lambda_b,x.lambda_behavior,x.step_psi,x.lambda2,inf_opts);
results = bayesopt(fun,[lambda_val,lambda_history,lambda_b,lambda_behavior,step_psi,lambda2],'Verbose',0,...
    'AcquisitionFunctionName','expected-improvement-plus','OutputFcn',@saveBOPlots);

%%
save(sprintf('bdlds_prmsrch_BayesOpt_%s.mat',todaysDate),'results')
%%
% https://www.mathworks.com/help/stats/bayesian-optimization-plot-functions.html
% https://www.mathworks.com/matlabcentral/answers/182574-save-all-the-plots
% https://www.mathworks.com/help/stats/bayesian-optimization-output-functions.html

function stop = saveBOPlots(results,state)
    stop = false;
    if strcmp(state, 'done')
        FolderName = 'bdLDSx_BOFigs';   % Your destination folder
        FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
        for iFig = 1:length(FigList)
          FigHandle = FigList(iFig);
          FigName   = get(FigHandle, 'Name');
          set(0, 'CurrentFigure', FigHandle);
          savefig(FigHandle, fullfile(FolderName, FigName, '.fig'));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%