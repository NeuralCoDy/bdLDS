%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

% parpool(16)

clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code
addpath(genpath('.'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
% rnndata   = load("./TaskRNNData/delayanti.mat"); 
% Test data (delayanti) has 80 trials
% 128 units (neurons)
% 125 time steps
% 20 inputs
% 3 outputs
rnnactivations = load("./TaskRNNData/share_mats/softplus/0/test/state.mat");
rnnoutput      = load("./TaskRNNData/share_mats/softplus/0/test/trial_output.mat");
rnninput       = load("./TaskRNNData/share_mats/softplus/0/test/trial_input.mat");
%%
todaysDate = input('Date (e.g., 251126):','s');

%% create cell arrays of neural and behavioral data
% set first dimension as units (neurons) for dFF, behavior features for
% behavior data
dFF          = num2cell(rnnactivations.delayanti, [1 3]);
dFF          = cellfun(@squeeze, dFF, 'UniformOutput', false);
dFF          = cellfun(@transpose, dFF, 'UniformOutput', false);
concatIO     = cat(3, rnninput.delayanti, rnnoutput.delayanti);
behaviorData = num2cell(concatIO, [1 3]);
behaviorData = cellfun(@squeeze, behaviorData, 'UniformOutput', false);
behaviorData = cellfun(@transpose, behaviorData, 'UniformOutput', false);
%% check data
% figure()
% tiledlayout('vertical')
% nexttile()
% plot(dFF{1})
% nexttile()
% plot(behaviorData{1}.')

%% check latent dimension (heuristic for setting nF here to start)
% 
% [~,~,~,~,explained,~] = pca(dFF{1}); % https://www.mathworks.com/matlabcentral/answers/713038-variance-explained-pca
% % rows: observations (time points), cols: variables (neurons)
% 
% figure();
% hold on
% bar(explained)
% plot(1:numel(explained), cumsum(explained), 'o-', 'MarkerFaceColor', 'r')
% yyaxis right
% h = gca;
% h.YAxis(1).Limits = [0 100];
% h.YAxis(1).Color = h.YAxis(1).Color;
% h.YAxis(1).TickLabel = strcat(h.YAxis(1).TickLabel, '%');
% h.YAxis(2).Limits = [0 100];
% h.YAxis(2).Color = h.YAxis(1).Color;
% h.YAxis(2).TickLabel = strcat(h.YAxis(2).TickLabel, '%');
% set(gcf,'Color','w')
% % title('Emmanuel - PCA variance explained, dF/F, no med filt or soft norm')
% title('Task RNN PCA variance explained')
% 
% hold off


% 8 dims --> 99% variance explained
% Double that
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Set parameters
inf_opts.nF              = 16; % # dynamics operators
inf_opts.N               = size(dFF{1},1); % # latent states
inf_opts.M               = size(dFF{1},1) ; % original dimension (# channels)
% disp('If reset dFF below, make sure to change inf_opts.M');
% inf_opts.lambda_val      = 0.9772;         % ASC added 7/25 % VARY - tradeoff: approx. SNR - see ratios of traces % 0.0001 (don't shrink this further - gets to be essentially 0 for solver_L1RLS)
% inf_opts.lambda_history  = 0.8702; % set from 500 iters BO %0.1;%0.01;         % ASC added 7/25
% inf_opts.lambda_b        = 0.2908;%0.25; %0.5;%0.39;           % ASC added 7/25 %0.025 % 0.5 %0.6
inf_opts.lambda_historyb = 1e-5;%0.45;%1.16;           % ASC added 7/25 %0.7 %0.45 %0.35 % 1.99 way too strong - all coefficients end up same
inf_opts.tol             = 1e-8;           % 1e-3
inf_opts.max_iter2       = 20;             % 500 %20
inf_opts.max_iters       = 10; %input('max iters:');
fprintf('\nmax_iters: %d\n',inf_opts.max_iters)
inf_opts.F_update        = true; % default = true;
inf_opts.D_update        = false; % NoObs case - default = true;
inf_opts.N_ex            = 1; %input('N_ex (e.g., 50):');%50;   
inf_opts.T_s             = size(dFF{1},2); %input('T_s (e.g., 200):');%200; 
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
inf_opts.lambda_f        = 1e-5;%0.1; %input('lambda_f (e.g., 0.2):'); %0.2 almost there but two operators nearly identical %0.4 too high - only 1-2 kinds of operators emerge
inf_opts.solver_type     = 'fista'; %input('Solver (fista or tfocs)?','s');%'tfocs'; % fista
inf_opts.special         = 'noobs';
inf_opts.deltaDynamics   = false; %default: false %use x_t-x_{t-1}
inf_opts.AcrossIndividuals     = 0; %input('Across individuals? 1 for yes, 0 for no:');

inf_opts.behaviordLDS    = 1; %input('behavior_dLDS?:');
if inf_opts.behaviordLDS == 1
    %inf_opts.step_psi        = 10; %input('step_psi:');%1; 
    inf_opts.lambda_behavior = 0.9313; % set from 500 iters BO %input('lambda_behavior:'); %1.6; %0.5521; % 0.5536 for opt for psi
    inf_opts.verysparsebhv   = 1;
    % inf_opts.lambda2         = 5; %input('lambda2:');
    inf_opts.psinorm         = 'frob'; %input('Norm on Psi - norm or frob?:','s');
else
    inf_opts.step_psi        = 0; % 30 for opt for psi
    inf_opts.lambda_behavior = 0; %0.5521; % 0.5536 for opt for psi
    inf_opts.verysparsebhv   = 0;
end


% inf_opts.solver_type     = 'fista'; %default: ''; (fista)
% inf_opts.CVX_Precision   = 'default';   %default = 'default'
% inf_opts.special         = ''; % regular bilinear inference
% inf_opts.debias          = true; % default = true; not used for NoObs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
rng('shuffle')
currentState01 = rng;
%% bayesopt

disp('If your results are too sparse (e.g., max b 0, median b use 0), set the upper bounds on the lambdas to smaller values')
disp('Something seems to be wrong on the server, run locally')
rng default
step_psi = optimizableVariable('step_psi',[1 20],'Type','integer'); 
lambda2 = optimizableVariable('lambda2',[0.1 10],'Type','real'); 
lambda_val = optimizableVariable('lambda_val',[0.01 0.3],'Type','real'); 
lambda_b = optimizableVariable('lambda_b',[0.01 0.3],'Type','real'); 
lambda_history = optimizableVariable('lambda_history',[0.01 0.3],'Type','real'); 
fun = @(x)runTaskRNNBO_opt5params_bdldsx(dFF,behaviorData,...
    x.step_psi,x.lambda2,x.lambda_val,x.lambda_b,x.lambda_history,inf_opts);
results = bayesopt(fun,[step_psi,lambda2,lambda_val,lambda_b,lambda_history],'Verbose',0,...
    'AcquisitionFunctionName','expected-improvement-plus','OutputFcn',@saveBOPlots);

%%
save(sprintf('taskRNN_prmsrch_BayesOpt_opt5params_bdldsx_%s.mat',todaysDate),'results')
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%