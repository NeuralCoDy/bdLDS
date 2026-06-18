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
    % Add all the required packages that are custom code
    addpath(genpath('../../npy-matlab/npy-matlab/'))                % Use steinmetz code for loady npy into matlab
    
    % addpath(genpath('../../../../home/adamsc/GITrepos/zebrafishDynamics/code/'))            % This is the main codebase for the project
    % addpath(genpath('../../../../home/adamsc/GITrepos/dynamics_learning/'))                 % This is the package for learning dynamics 
    
    addpath(genpath('../../../adamsc/data/zebrafish/'))

    dat.dir   = './FishBehaviorData/';
    dat.Ffile = 'for_JHU_cells_dff.npy';
    dat.Bfile = 'for_JHU_gad1b-6hz-backwardpulse_integrate_behavior_ds.npy';
    dat.Vfile = 'for_JHU_-gad1b-6hz-backwardpulse_integrate_visual_velocity_ds.npy';
else
    addpath(genpath('.'))
    load("saveFish_En_250901_lbhv1_5000i.mat");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data

if ifForCIS
    dFF = readNPY([dat.dir, dat.Ffile]);   % getNPYarray doesn't exist % Neuron (rows) x time matrix of fluorescent activity
    bhv = readNPY([dat.dir, dat.Bfile]);                               % 1 x time vector of behavioral data
    vel = readNPY([dat.dir, dat.Vfile]);                               % 1 x time vector of velocity
end
%%
dFFOriginal = dFF;

[depVar,dFF,timePoints] = preprocessABehaviorVariable(dFF,bhv,vel,'bhv','later');
% [depVar2,~,~]            = preprocessABehaviorVariable(dFF,bhv,vel,'vel','later');

ifDoBoutType = input('Include bout type information in behavior? (1 or 0):');

if ifDoBoutType
    [depVar3,dFF,timePoints] = preprocessABehaviorVariable(dFFOriginal,bhv,vel,'boutType','later');
    warning('If use boutType, does not use last incomplete bout, so 4200 timepoints instead of 4282')
end

% depVar  = horzcat(depVar1,depVar2,depVar3);
% depVar  = horzcat(depVar1,depVar2);

doSmoothedBehavior = input('Smooth behavior?');

if doSmoothedBehavior
    smoothingwindow = input('Sigma for smoothing?');
    sigma           = smoothingwindow; 
    gaussian_range  = -3*sigma:3*sigma; % setting up Gaussian window
    gaussian_kernel = normpdf(gaussian_range,0,sigma); % setting up Gaussian kernel
    gaussian_kernel = gaussian_kernel/sum(gaussian_kernel);
    depVar          = conv(depVar(:),gaussian_kernel(:),'same');
    
end

depVar = depVar./(max(abs(depVar))); %1D behavior

if ifDoBoutType
    depVar = horzcat(depVar(1:4200,:),depVar3);
end

disp(size(depVar))

%% deal with NaNs
[row, col] = find(isnan(dFF));
dFF(row,:) = [];
% dFF(isnan(dFF)) = 0;

%%
if ~ifForCIS
    figure()
    set(gcf,'Color','w')
    histogram(dFF)
    title('Neuron data - unprocessed')
    
    figure()
    set(gcf,'Color','w')
    imagesc(dFF)
    colorbar
    xlabel('Neurons')
    xlabel('Time')
    title('Neuron data - unprocessed')
    
    figure()
    set(gcf,'Color','w')
    plot(bhv)
    xlabel('Time')
    title('Behavior')
end

%% select and set dFF 
dFFOriginal = dFF;
%%
[filtered.dFF] = filterByStd(dFFOriginal,1,0,0);
%%
if ~ifForCIS
    figure()
    tiledlayout("flow")
    nexttile
    set(gcf,'Color','w')
    histogram(filtered.dFF)
    title('En data, capped')

    stdFilt = robustSTD(filtered.dFF.');
    nexttile
    set(gcf,'Color','w')
    histogram(stdFilt)
    title('En stdev, capped')
    
    nexttile
    set(gcf,'Color','w')
    imagesc(filtered.dFF)
    xlabel('Neurons')
    xlabel('time points')
    title('En data, capped')
    colorbar
end
%%
dFF = filtered.dFF;

%titleForPCAFig = "Filtered dFF 5 std, By Neuron, Zero, Mode normalized";

% trialIdentifier = input("Trial identifier:", "s");

%% check - should be no NaNs
if sum(isnan(dFF),'all') > 0
    warning('NaNs remain - check preprocessing');
end
%% check - should be no Infs
if sum(isinf(dFF),'all') > 0
    warning('Infs remain - check preprocessing');
end
%% check - each row should have a non-zero std (otherwise capping was too aggressive)
if sum(std(dFF,0,2)==0,'all') > 0
    warning('Cap too aggressive, data all one value in some rows - check preprocessing');
end
%% median filtering

doMedianFilter = input('Median filter (smoothing)? (1 or 0):');
if doMedianFilter
    medFiltSize = input('Median filter window size (recommended > 5):');
    dFF = medfilt2(dFF,[1,medFiltSize]);
end

%% rescale data to mode (by neuron) after filtering

doDFOverF = input('Do dF/F (1) or just dF (0)?:');

modeEst      = kernelModeEstimate(dFF.');
dFFMinusMode = bsxfun(@minus, dFF.', modeEst);

if doDFOverF
    stdDFF       = robustSTD(dFF.');
    softnormterm = input('Soft norm term (rec. 1):'); %prctile(stdDFF,50);
    dFFModeNorm  = bsxfun(@times, dFFMinusMode, 1./(stdDFF+softnormterm)); 
    disp('Soft normalization on std')
    dFF          = dFFModeNorm.';
else
    dFF          = dFFMinusMode.';
end



% dFFModeRescaled = dFFModeNorm.';

%%
if ~ifForCIS
    if doDFOverF

        disp('update titles')

        figure()
        tiledlayout("flow")
        nexttile
        set(gcf,'Color','w')
        histogram(dFF)
        % title('En data, no med filt, no soft norm, dFF')
        title('En data, med filt 5, soft norm +1, dFF')
    
        postDFFstd = robustSTD(dFF.');
        nexttile
        set(gcf,'Color','w')
        histogram(postDFFstd)
        % title('En stdev, no med filt, no soft norm, dFF')
        title('En stdev, med filt 5, soft norm +1, dFF')

        
        nexttile
        set(gcf,'Color','w')
        imagesc(dFF)
        xlabel('Neurons')
        xlabel('time points')
        % title('En data, no med filt, no soft norm, dFF')
        title('En data, med filt 5, soft norm +1, dFF')
        colorbar
    end
end
% figure()
% set(gcf,'Color','w')
% histogram(dFF)
% title('Neuron data - preprocessed, rescaled')
% 
% figure()
% set(gcf,'Color','w')
% imagesc(dFF)
% xlabel('Neurons')
% xlabel('Time')
% title('Neuron data - preprocessed, rescaled')
% colorbar

%% create cell array with your samples 

%Examples

% % 1 worm, 1 trial
% whichStimWorm = 1;
% dFF = dFF_S_AllWorms(whichStimWorm);
% inf_opts.AcrossIndividuals     = false;

% % 1 worm, multiple durations of trials (artificially trimmed) (same channels, could be different trial durations)
% whichStimWorm = 1;
% dFF = {dFF_S_AllWorms{whichStimWorm}(:,1:2000);dFF_S_AllWorms{whichStimWorm}(:,1:1000);dFF_S_AllWorms{whichStimWorm}(:,1:500)}; 
% bhv = {bhv_S_AllWorms{whichStimWorm}(:,1:2000);bhv_S_AllWorms{whichStimWorm}(:,1:1000);bhv_S_AllWorms{whichStimWorm}(:,1:500)}; 

% % multiple worms (different channels, different trial durations)
% dFF = dFF_S_AllWorms;
% bhv = bhv_S_AllWorms;
% inf_opts.AcrossIndividuals     = true;

% 1 fish, 1 trial
% inf_opts.AcrossIndividuals     = false;
% dFF{1}

% behavior dLDS
% inf_opts.behaviordLDS = true;
%inf_opts.lambda_behavior = 0.1;
%inf_opts.nBhv = size(behavior_data,2);



%%
% for ii = 1:size(dFF,1)
%     thisDFF = dFF{ii};
%     thisDFF(thisDFF<0) = 0;
%     thisDFF = thisDFF./max(thisDFF,[],'all');
%     dFF{ii} = thisDFF;
% end

%% test PCA dims

if ~ifForCIS
    disp('update titles')

    % [~,~,~,~,explained,~] = pca(dFFOriginal.'); % https://www.mathworks.com/matlabcentral/answers/713038-variance-explained-pca
    % % rows: observations (time points), cols: variables (neurons)
    % 
    % figure();
    % hold on
    % bar(explained)
    % plot(1:numel(explained), cumsum(explained), 'o-', 'MarkerFaceColor', 'r')
    % yyaxis right
    % h = gca;
    % h.YAxis(2).Limits = [0 100];
    % h.YAxis(2).Color = h.YAxis(1).Color;
    % h.YAxis(2).TickLabel = strcat(h.YAxis(2).TickLabel, '%');
    % set(gcf,'Color','w')
    % title('PCA variance explained')
    % hold off
    
    [~,~,~,~,explained,~] = pca(dFF); % https://www.mathworks.com/matlabcentral/answers/713038-variance-explained-pca
    % rows: observations (time points), cols: variables (neurons)
    
    figure();
    hold on
    bar(explained)
    plot(1:numel(explained), cumsum(explained), 'o-', 'MarkerFaceColor', 'r')
    yyaxis right
    h = gca;
    h.YAxis(1).Limits = [0 100];
    h.YAxis(1).Color = h.YAxis(1).Color;
    h.YAxis(1).TickLabel = strcat(h.YAxis(1).TickLabel, '%');
    h.YAxis(2).Limits = [0 100];
    h.YAxis(2).Color = h.YAxis(1).Color;
    h.YAxis(2).TickLabel = strcat(h.YAxis(2).TickLabel, '%');
    set(gcf,'Color','w')
    title('En - PCA variance explained, dF/F, no med filt or soft norm')
    hold off
    
    % pcadim = 2;
    % while sum(explained(1:pcadim)) < 0.95
    %     pcadim = pcadim + 1;
    % end
    % 
    % disp(pcadim-1);
    % disp(sum(explained(1:pcadim-1)))
end

%%
% dFF =  dFFMinusMode.'; %filtered.dFF;
% dFF = dFFModeNorm.';

%%
doLB = input('With Ljung-Box? (1 or 0):');

if doLB
    lbPerNeuronOriginal = zeros(1,size(dFF,1)); % here, dFF is not a cell array
    
    singleDFF = dFF; % save RAM
    
    a = double(singleDFF);
    
    parfor jj = 1:size(a,1)
        [lbPerNeuronOriginal(1,jj),...
            ~,...
            ~,~]...
            = lbqtest1d(a(jj,:),'Lags',[500 1000]);
     
    end
    
    dFF = dFF(lbPerNeuronOriginal(1,:)==1,:);
end
%%
doRescaleToMaxAbs = input('Rescale to max abs value? (0 or 1):');
if doRescaleToMaxAbs
    all_max = max(abs(dFF(:)));
    dFF = dFF.*1./all_max;
end
% all_max = cell2mat(cellfun(@(x)max(abs(x(:))),dFF,'UniformOutput',false));
% all_max_max = max(all_max(:));
% 
% for ii = 1:size(dFF,1)
%     thisDFF = dFF{ii};
%     % thisDFF(thisDFF<0) = 0; - not relevant - no negative rates here
%     thisDFF = thisDFF.*1./all_max_max;
%     dFF{ii} = thisDFF;
% end
%%
if ~ifForCIS
    [~,~,~,~,explained,~] = pca(dFF); % https://www.mathworks.com/matlabcentral/answers/713038-variance-explained-pca
    % rows: observations (time points), cols: variables (neurons)
    
    figure();
    hold on
    bar(explained)
    plot(1:numel(explained), cumsum(explained), 'o-', 'MarkerFaceColor', 'r')
    yyaxis right
    h = gca;
    h.YAxis(2).Limits = [0 100];
    h.YAxis(2).Color = h.YAxis(1).Color;
    h.YAxis(2).TickLabel = strcat(h.YAxis(2).TickLabel, '%');
    set(gcf,'Color','w')
    title('PCA variance explained')
    hold off

    figure(); 
    tiledlayout('flow')
    nexttile()
    histogram(dFF(:)); title('Data')
    
    disp(var(dFF,0,2))
    
    nexttile(); histogram(var(dFF,0,2)); title('Neuron var')
    
    fanoFactorish = var(dFF,0,2)/mean(dFF,2);
    nexttile();histogram(fanoFactorish);title('Neuron var/mean')
    
    disp(prctile(dFF,[5 95],2))
    
    nexttile(); histogram(prctile(dFF,5,2)); title('5th percentile')
    nexttile(); histogram(prctile(dFF,95,2)); title('95th percentile')

end
%% FIXME: add in downsampling to 1 Hz from 6 Hz
% every6 = 1:6:size(dFF,2);
% dFF = dFF(:,every6);
%%
% fftDFF = fft(dFFMinusMode);
% fftDFFMag = abs(fftshift(fftDFF));
% %%
% figure()
% histogram(fftDFFMag)
%%
if ~iscell(dFF)
    dFFcell{1} = dFF;
    dFF = dFFcell;
end

checkBehaviorDims= input('How many feature dimensions should your behavior have?');

if size(depVar,2) == checkBehaviorDims
    behavior_data = depVar.'; %N by T
elseif size(depVar,1) == checkBehaviorDims
    disp('No transpose needed: rows are features, columns are time points')
    behavior_data = depVar;
else
    warning('Check your behavior dimensions')
end
    

if ~iscell(behavior_data)
    bhvcell{1} = behavior_data;
    behavior_data = bhvcell;
end

inf_opts.AcrossIndividuals     = false;

if inf_opts.AcrossIndividuals
    [inf_opts.M,~]               = cellfun(@size,dFF,'UniformOutput',false);%size(dFF,1) ; % EY changed 9/6/23 for ManyWorms case
else
    inf_opts.M = size(dFF{1},1) ; % original dimension (# channels), N by T
end

%%
% rng default
% lambda_behavior = optimizableVariable('lambda_behavior',[0.01 2],'Type','real'); % [0.001 1] --> 0.5536
% lambda_b = optimizableVariable('lambda_b',[0.01 2],'Type','real'); % [0.1 10], real --> 9.9994
% fun = @(x)runBehaviorDLDSFish_LambdaBehaviorLambdaB_BestYhat(dFF,behavior_data,...
%     x.lambda_behavior,x.lambda_b);
% results = bayesopt(fun,[lambda_behavior,lambda_b],'Verbose',0,...
%     'AcquisitionFunctionName','expected-improvement-plus');

%%
% Set parameters
inf_opts.nF              = input('nF (sugg. > 10):'); % 10 not enough - using all
inf_opts.N               = input('n:'); %scaled down from 100 (zebrafish, whole brain) 

disp('If reset dFF below, make sure to change inf_opts.M');
inf_opts.lambda_val      = input('lambda_val (sugg. 0.1):'); % 0.1 % VARY - tradeoff: approx. SNR - see ratios of traces
inf_opts.lambda_history  = input('lambda_history (sugg. 0.9):'); % 0.9
inf_opts.lambda_b        = input('lambda_b (sugg. 0.01):'); % 0.01 %FIXME: stronger than 0.01 (also other reg params) - 0.1 all disappear, 0.05 only 2 traces
inf_opts.lambda_historyb = 0;
inf_opts.behaviordLDS    = input('behavior? (0 or 1):');
inf_opts.lambda_behavior = input('lambda_behavior (sugg. 1.3):');
inf_opts.tol             = 1e-3; % 1e-3
inf_opts.max_iter2       = 50; %500 
inf_opts.max_iters       = input('max_iters (sugg. > 2000):');
inf_opts.special         = '';
inf_opts.F_update        = true;
inf_opts.D_update        = true;
% inf_opts.AcrossIndividuals     = false;
% inf_opts.AcrossIndividuinals     = 0; %input('Across individuals? 1 for yes, 0 for no:');
inf_opts.N_ex            = input('N_ex:');%40;
inf_opts.T_s             = input('T_s:');%30;
% inf_opts.sampleProportionally = true; %default: prop by trial length
% inf_opts.special         = 'noobs';
inf_opts.solver_type          = input('Solver (tfocs or fista):','s');
inf_opts.nBhv            = input('nBhv:');



% inf_opts.step_d          = 20;  % ASC added 7/25
% inf_opts.step_f          = 30;  % ASC added 7/25
% inf_opts.plot_option     = 20;  % ASC added 7/25
inf_opts.lambda_f        = input('lambda_f:'); % ASC added 7/28'
inf_opts.step_psi        = 10;

% inf_opts.verysparsebhv   = input('very sparse suspected? (1 or 0):');
% 

%%
% if ifForCIS
%     addpath(genpath('~/my_documents/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code'))
% end

currentState = rng;

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


varExpl = overallVarExpl(dFF,A_cell,Phi,inf_opts);
%%
varExplPerNeuron = perNeuronVarExpl(dFF,A_cell,Phi,inf_opts);
%% look at model output traces - coefficients, reconstruction
trialIdentifier = input('Trial identifier:','s');
dataIdentifier  = 'FishBehavior';
samplingRate    = 6;
%%
plotSummaryCandReconstr_OneC(dFF,Phi,B_cell,A_cell,inf_opts,trialIdentifier,varExplPerNeuron,samplingRate,dataIdentifier);
% plotSummaryCandReconstr_OneC(dFF,Phi,B_cell,A_cell,inf_opts,trialIdentifier,varExplPerNeuron,samplingRate,dataIdentifier,[],[],[1 3 7 15 22 23]); % [3 4 22 23]
%%
plotSummaryCandReconstr(dFF,Phi,B_cell,A_cell,inf_opts,trialIdentifier,varExplPerNeuron,samplingRate,dataIdentifier);
%%
ifDisplayOnSecondMonitor = input('Display on second monitor?(1 yes, 0 no):');
p = get(0, "MonitorPositions");
c_values = B_cell{1}.';
Phit = Phi.';

for timePoint = 200 %[200 400 600 800 1000 1200 1400]
    f = figure();
    if ifDisplayOnSecondMonitor
        f.Position = p(2,:);
    else
        f.Position = p(1,:);
    end
    f.WindowState = 'maximized';
    
    scaleEachDOByC = F;
    for whichDO = 1:size(F,2)
        scaleEachDOByC{whichDO} = c_values(timePoint,whichDO).*F{whichDO};
    end
    scaleEachDOByCArray = cat(3,scaleEachDOByC{:});
    dynsTogether =  sum(scaleEachDOByCArray,3);
    [U,L]   = eig(dynsTogether);
    maxEvalScaleFactor = max(abs(L(:)));
    scaledF = real(U*L*(diag(1./maxEvalScaleFactor)*inv(U)));
    
    hold on;
    connectivityMatrix                          = Phi * scaledF * Phit;
    connThresholded                             = connectivityMatrix;
    upperQuantileThresholdStrength              = quantile(abs(connectivityMatrix),0.9999,"all");
    connThr                                     = upperQuantileThresholdStrength; %0.18 * max(abs(connectivityMatrix(:))); % 1 stdev
    connThresholded(abs(connThresholded) < connThr) = 0;
    connThresholdedAbs = abs(connThresholded);
    circularGraph(connThresholdedAbs);
    % circularGraph(connThresholdedAbs,...
    %     'Label',cellstr(channelNamesStrings));

    hold off
end
%%

find(abs(Psi) > 0.2)
find(abs(Psi) < 0.05)
%%
figure()
for i = 1:25 %[8 10 18 25 2 5 16 24]%1:size(B_cell{1})
    for j = 1:inf_opts.nBhv
        cla
    
        nTimepoints = size(A_cell{1},2);
        tVals = linspace(0,nTimepoints/samplingRate,nTimepoints);
    
        % tiledlayout(4,1)
        tiledlayout('vertical')
    
        nexttile
        plot(tVals,B_cell{1}(i,:))
        box off
        axis([0 max(tVals) -Inf Inf]);
        rBhvCo = corrcoef(B_cell{1}(i,:),behavior_data{1}(j,:));
        title(sprintf('Dyn coef %d, Psi coef %0.2f, R2 vs behavior %d %0.2f',i,Psi(j,i), j, rBhvCo(1,2)^2))
    
        nexttile
        % plot(behavior_data{1}(2,:))
        plot(tVals,behavior_data{1})
        box off
        axis([0 max(tVals) -Inf Inf]);
        title('Behavior (smoothed), input to model')
    
    
        nexttile
        bhvReconstruction = Psi*B_cell{1};
        plot(tVals,bhvReconstruction)
        box off
        axis([0 max(tVals) -Inf Inf]);
        title('Behavior reconstruction')
        legend(['motor (1)'; '2        '; '3        '; '4        '])
        xlabel('Time (s)')
    
    
        if doSmoothedBehavior
            sigma           = smoothingwindow; 
            gaussian_range  = -3*sigma:3*sigma; % setting up Gaussian window
            gaussian_kernel = normpdf(gaussian_range,0,sigma); % setting up Gaussian kernel
            gaussian_kernel = gaussian_kernel/sum(gaussian_kernel);
            depVarRe        = conv(bhvReconstruction(:),gaussian_kernel(:),'same');
    
        end
    
        % depVarRe = depVarRe./(max(abs(depVarRe))); %1D behavior
    
        % nexttile
        % plot(tVals,depVarRe)
        % box off
        % axis([0 max(tVals) -Inf Inf]);
        % title('Behavior reconstruction (smoothed)')
    
    
        rBhvRe = corrcoef(depVarRe,behavior_data{1});
        % rBhvCo = corrcoef(B_cell{1}(i,:),behavior_data{1}(2,:));
    
        disp(i)
        disp(rBhvRe(1,2)^2)
        disp(rBhvCo(1,2)^2)
        disp(Psi(i))
    
        pause()

    end
end
%%
whichData = varExplPerNeuron;
histVarExpl(whichData,inf_opts,trialIdentifier,dataIdentifier);
%%
mapFishHindbrainData(varExplPerNeuron,lbPerNeuronOriginal);
%%
mapFishHindbrainData(dFF{1},lbPerNeuronOriginal);
%%
% traceVsBlock = lookupTraceVsBlock(lbPerNeuronOriginal,64,70);
% traceVsBlock = lookupTraceVsBlock(lbPerNeuronOriginal,1,1);
traceVsCentroid = lookupTraceVsCentroid(lbPerNeuronOriginal);

c_values = B_cell{1}.';
%%
% mapFishHindbrainConnectivity(varExplPerNeuron,lbPerNeuronOriginal,Phi,F,c_values,0,traceVsBlock,64,70)
% mapFishHindbrainConnectivity(varExplPerNeuron,lbPerNeuronOriginal,Phi,F,c_values,0,traceVsBlock,1,1)
%% operators
mapFishHindbrainConnectivity(varExplPerNeuron,lbPerNeuronOriginal,Phi,F,c_values,0,traceVsCentroid,trialIdentifier); %do operators
%% time points
mapFishHindbrainConnectivity(varExplPerNeuron,lbPerNeuronOriginal,Phi,F,c_values,1,traceVsCentroid,trialIdentifier); %do timepoints

%%
if inf_opts.AcrossIndividuals
    for ii = 1:size(Phi,1)
        x_values = cell2mat(A_cell{ii});
        dataReconstructed = (Phi{ii} * x_values).';    
        
        flatDFF = reshape(dFF{ii}, numel(dFF{ii}),[]);
        flatReconstructed = reshape(dataReconstructed,numel(dataReconstructed),[]);
        rval = corrcoef(dFF{ii}, dataReconstructed.');
        varExpl = (rval(1,2))^2;
        disp(varExpl);

    end
else
    for ii = 1:size(A_cell,1) %in the event of multiple trials
        x_values = A_cell{ii};
        dataReconstructed = (Phi * x_values).';    
        
        maxvalDFF = max(dFF{ii},[],'all');
        stdzdOn = false;
        if stdzdOn
            dataReconstructedRescaled = maxvalDFF .*(Phi * x_values);
            dFFRescaled = maxvalDFF .* dFF{ii};
        else
            dFFRescaled = dFF{ii};
            dataReconstructedRescaled = (Phi * x_values);
        end
       
        
        flatDFF = reshape(dFFRescaled, numel(dFFRescaled),[]);
        flatReconstructed = reshape(dataReconstructedRescaled,numel(dataReconstructedRescaled),[]);
        rval = corrcoef(flatDFF, flatReconstructed);
        varExpl = (rval(1,2))^2;
        disp(varExpl);
        % pause();
    end
end
%%
[~,whichTracesMin] = mink(varExplPerNeuron(:,ii),1000);
[~,whichTracesMax] = maxk(varExplPerNeuron(:,ii),1000);

%%
figure()
tiledlayout(4,1)
nexttile
imagesc(dFF{1}(whichTracesMin,:))
title('orig worst')
colorbar
nexttile
imagesc(dataReconstructedRescaled(whichTracesMin,:))
colorbar
title('recon worst')

nexttile
imagesc(dFF{1}(whichTracesMax,:))
colorbar
title('orig best')
nexttile
imagesc(dataReconstructedRescaled(whichTracesMax,:))
colorbar
title('recon best')
%%
% Fs = 66;            % Sampling frequency                    
% T  = 1/Fs;             % Sampling period       
% L  = 4282;             % Length of signal
% t  = (0:L-1)*T;        % Time vector
% 
% %note: FFT goes by columns
% 
% figure()
% tiledlayout(4,1)
% ax1 = nexttile;
% X = dFF{1}(whichTracesMin,:);
% Y = fft(X.');
% % plot(Fs/L*(0:L-1),abs(Y),"LineWidth",3)
% % title("Complex Magnitude of fft Spectrum - orig data worst recon")
% % xlabel("f (Hz)")
% % ylabel("|fft(X)|")
% plot(Fs/L*(-L/2:L/2-1),abs(fftshift(Y)),"LineWidth",3)
% title("fft Spectrum in the Positive and Negative Frequencies- orig worst")
% xlabel("f (Hz)")
% ylabel("|fft(X)|")
% ax2 = nexttile;
% X = dataReconstructedRescaled(whichTracesMin,:);
% Y = fft(X.');
% % plot(Fs/L*(0:L-1),abs(Y),"LineWidth",3)
% % title("Complex Magnitude of fft Spectrum - worst recon")
% % xlabel("f (Hz)")
% % ylabel("|fft(X)|")
% plot(Fs/L*(-L/2:L/2-1),abs(fftshift(Y)),"LineWidth",3)
% title("fft Spectrum in the Positive and Negative Frequencies- recon worst")
% xlabel("f (Hz)")
% ylabel("|fft(X)|")
% 
% ax3 = nexttile;
% X = dFF{1}(whichTracesMax,:);
% Y = fft(X.');
% % plot(Fs/L*(0:L-1),abs(Y),"LineWidth",3)
% % title("Complex Magnitude of fft Spectrum - orig data best recon")
% % xlabel("f (Hz)")
% % ylabel("|fft(X)|")
% plot(Fs/L*(-L/2:L/2-1),abs(fftshift(Y)),"LineWidth",3)
% title("fft Spectrum in the Positive and Negative Frequencies- orig best")
% xlabel("f (Hz)")
% ylabel("|fft(X)|")
% ax4 = nexttile;
% X = dataReconstructedRescaled(whichTracesMax,:);
% Y = fft(X.');
% % plot(Fs/L*(0:L-1),abs(Y),"LineWidth",3)
% % title("Complex Magnitude of fft Spectrum - best recon")
% % xlabel("f (Hz)")
% % ylabel("|fft(X)|")
% plot(Fs/L*(-L/2:L/2-1),abs(fftshift(Y)),"LineWidth",3)
% title("fft Spectrum in the Positive and Negative Frequencies- recon best")
% xlabel("f (Hz)")
% ylabel("|fft(X)|")
% linkaxes([ax1 ax2 ax3 ax4])
% 
% 
% %%
% figure()
% orig = dFF{1}(whichTracesMin,:);
% recon = dataReconstructedRescaled(whichTracesMin,:);
% tiledlayout(4,1)
% ax1 = nexttile;
% histogram(orig(:))
% title('orig worst')
% ax2 = nexttile;
% histogram(recon(:))
% title('recon worst')
% 
% orig = dFF{1}(whichTracesMax,:);
% recon = dataReconstructedRescaled(whichTracesMax,:);
% ax3 = nexttile;
% histogram(orig(:))
% title('orig best')
% ax4 = nexttile;
% histogram(recon(:))
% title('recon best')
% linkaxes([ax1 ax2 ax3 ax4])
% 
% %%
% figure()
% orig = std(dFF{1}(whichTracesMin,:),[],2);
% recon = std(dataReconstructedRescaled(whichTracesMin,:),[],2);
% tiledlayout(4,1)
% ax1 = nexttile;
% histogram(orig)
% title('std orig worst')
% ax2 = nexttile;
% histogram(recon)
% title('std recon worst')
% 
% orig = std(dFF{1}(whichTracesMax,:),[],2);
% recon = std(dataReconstructedRescaled(whichTracesMax,:),[],2);
% ax3 = nexttile;
% histogram(orig)
% title('std orig best')
% ax4 = nexttile;
% histogram(recon)
% title('std recon best')
% linkaxes([ax1 ax2 ax3 ax4])
% 
% %%
% figure()
% 
% orig = dFF{1}(whichTracesMin,:);
% recon = dataReconstructedRescaled(whichTracesMin,:);
% tiledlayout(4,1)
% ax1 = nexttile;
% histogram(std(abs(fftshift(fft(orig.')))))
% title('orig worst')
% ax2 = nexttile;
% histogram(std(abs(fftshift(fft(recon.')))))
% title('recon worst')
% 
% orig = dFF{1}(whichTracesMax,:);
% recon = dataReconstructedRescaled(whichTracesMax,:);
% ax3 = nexttile;
% histogram(std(abs(fftshift(fft(orig.')))))
% title('orig best')
% ax4 = nexttile;
% histogram(std(abs(fftshift(fft(recon.')))))
% title('recon best')
% linkaxes([ax1 ax2 ax3 ax4])
% 
% %%
% figure
% for i = 1:size(whichTracesMin,1)
%     cla
%     tiledlayout("flow")
%     nexttile
%     autocorr(double(dFF{1}(whichTracesMin(i),:)),100)
%     [h,pValue,stat,cValue] = lbqtest(double(dFF{1}(whichTracesMin(i),:)),'Lags',100);
%     pValue
%     h
%     nexttile
%     autocorr(double(dataReconstructedRescaled(whichTracesMin(i),:)),100)
%     nexttile
%     autocorr(double(dFF{1}(whichTracesMax(i),:)),100)
%     nexttile
%     autocorr(double(dataReconstructedRescaled(whichTracesMax(i),:)),100)
% 
%     pause()
% end
% %%
% [lbPerNeuronOriginalSMin,lbPerNeuronReconstrSMin,lbPerNeuronMinusReconstrSMin] = perNeuronLjungBox(dFF,A_cell,Phi,inf_opts,whichTracesMin);
% %%
% [lbPerNeuronOriginalSMax,lbPerNeuronReconstrSMax,lbPerNeuronMinusReconstrSMax] = perNeuronLjungBox(dFF,A_cell,Phi,inf_opts,whichTracesMax);
% %% All means after all preprocessing
% [lbPerNeuronOriginalSAll,lbPerNeuronReconstrSAll,lbPerNeuronMinusReconstrSAll] = perNeuronLjungBox(dFF,A_cell,Phi,inf_opts);
% %% plot behavior vs. reconstruction
% dBhv = [0 diff(bhvcell{1})];
% reconstructedBehavior = Psi*B_cell{1};
% reconstructeddBhv = [0 diff(reconstructedBehavior)];
% rval = corrcoef(bhvcell{1}, reconstructedBehavior);
% varExplBhv = (rval(1,2))^2;
% rval = corrcoef(dBhv, reconstructeddBhv);
% varExpldBhv = (rval(1,2))^2;
% figure()
% tiledlayout(2,1)
% nexttile
% plot(bhvcell{1})
% nexttile
% plot(reconstructedBehavior)
% % title(sprintf("behavior R2 %.2f",varExplBhv))
% % nexttile
% % plot(dBhv)
% % nexttile
% % plot([0 diff(reconstructedBehavior)])
% % title(sprintf("dBhv R2 %.2f",varExpldBhv))
% 
% %%
% filename = input("workspace filename:","s");
% save(filename);
% 
% 
% 
% for ii = 1% 1:size(Phi,1)
%     fig1 = figure();
%     subplot(3,1,1)
%     c_values = cell2mat(B_cell{ii});
%     plot(c_values.')
%     hold on
%     ylabel("c")
%     hold off
% 
%     subplot(3,1,2)
%     x_values = cell2mat(A_cell{ii});
%     plot((Phi{ii} * x_values).')
%     hold on
%     ylabel("reconstr. fluor.")
%     hold off
% 
%     subplot(3,1,3)
%     plot(dFF{ii}.')
%     hold on
%     ylabel("scaled fluor.")
%     hold off
% end
% 
% filename = input("figure filename:","s");
% saveas(fig1,filename);
% close all
% 
% ids = dat.WT_Stim(whichStimWorm).IDs;
% identifiedNeurons = find(cellfun(@isempty,ids)==0);
% neurSel = randsample(identifiedNeurons, 10);

%%
function [lbPerNeuronOriginalS,lbPerNeuronReconstrS,lbPerNeuronMinusReconstrS] = perNeuronLjungBox(dFF,A_cell,Phi,inf_opts,varargin)
if nargin > 4
    whichIndices = varargin{1};
else
    whichIndices = 1:size(dFF{1},1);
end
if inf_opts.AcrossIndividuals
    disp('These dims can be different for every individual - must'+...+
        'be cell')
else
    lbPerNeuronOriginal        = zeros(size(dFF{1}(whichIndices,:),1),size(A_cell,1)); % neurons by trials
    lbPerNeuronReconstr        = zeros(size(dFF{1}(whichIndices,:),1),size(A_cell,1)); % neurons by trials
    lbPerNeuronMinusReconstr   = zeros(size(dFF{1}(whichIndices,:),1),size(A_cell,1)); % neurons by trials

    lbPerNeuronOriginalP       = zeros(size(dFF{1}(whichIndices,:),1),size(A_cell,1)); % neurons by trials
    lbPerNeuronReconstrP       = zeros(size(dFF{1}(whichIndices,:),1),size(A_cell,1)); % neurons by trials
    lbPerNeuronMinusReconstrP  = zeros(size(dFF{1}(whichIndices,:),1),size(A_cell,1)); % neurons by trials

    lbPerNeuronOriginalS       = zeros(size(dFF{1}(whichIndices,:),1),size(A_cell,1)); % neurons by trials
    lbPerNeuronReconstrS       = zeros(size(dFF{1}(whichIndices,:),1),size(A_cell,1)); % neurons by trials
    lbPerNeuronMinusReconstrS  = zeros(size(dFF{1}(whichIndices,:),1),size(A_cell,1)); % neurons by trials

    for ii = 1:size(A_cell,1) %in the event of multiple trials
        disp(ii)
        x_values = A_cell{ii};

        singleDFF = dFF{ii}; % save RAM
        singleReconstr = (Phi * x_values);

        a = double(singleDFF);
        b = double(singleReconstr);

        a = a(whichIndices,:);
        b = b(whichIndices,:);

        parfor jj = 1:size(a,1)
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

            [lbPerNeuronOriginal(jj,ii),...
                lbPerNeuronOriginalP(jj,ii),...
                lbPerNeuronOriginalS(jj,ii),~]...
                = lbqtest1d(a(jj,:),'Lags',[500 1000]);
            [lbPerNeuronReconstr(jj,ii),...
                lbPerNeuronReconstrP(jj,ii),...
                lbPerNeuronReconstrS(jj,ii),~]...
                = lbqtest1d(b(jj,:),'Lags',[500 1000]);
            [lbPerNeuronMinusReconstr(jj,ii),...
                lbPerNeuronMinusReconstrP(jj,ii),...
                lbPerNeuronMinusReconstrS(jj,ii),~]...
                = lbqtest1d(a(jj,:)-b(jj,:),'Lags',[500 1000]);

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
    end
end

end

%%
function [depVar, dFF, timePoints] = preprocessABehaviorVariable(dFF,bhv,vel,whichDepVar,whichPoints)

velTrimmed = vel(92:5041); % https://www.mathworks.com/matlabcentral/answers/326675-cumulative-sum-at-over-an-specified-interval
B=reshape(velTrimmed,150,[]);
C=cumsum(B,1);
displacementPerBout=reshape(C,[],1);

% whichDepVar = input('Which dependent variable? bhv,cumsumbhv,vel,cumsumvel,displBout,boutType:','s');
% whichPoints = input('Later points only, or all points? later, all:','s');
% 
% dataIdentifier = "Zebrafish";
% trialIdentifier = input('Supply a unique trial identifier string:','s');
% if isempty(trialIdentifier)
%     trialIdentifier = string(randi(100000000));
% end

pooledPrediction = input('Pooled:');
if isempty(pooledPrediction)
    pooledPrediction = 0;
end

stackedPrediction = input('Stacked:');
if isempty(stackedPrediction)
    stackedPrediction = 0;
end

timePoints = 1:size(dFF,2);

if strcmp(whichDepVar,'bhv')
    depVar = bhv;
    if pooledPrediction|stackedPrediction
        depVar = depVar(92:5041);
        dFF = dFF(:,92:5041); % trimmed for fair comparison
        timePoints = 92:5041;
    end
elseif strcmp(whichDepVar,'cumsumbhv')
    depVar = cumsum(bhv);
    if pooledPrediction|stackedPrediction
        depVar = depVar(92:5041);
        dFF = dFF(:,92:5041); % trimmed for fair comparison
        timePoints = 92:5041;
    end
elseif strcmp(whichDepVar,'vel')
    depVar = vel;
    if pooledPrediction|stackedPrediction
        depVar = depVar(92:5041);
        dFF = dFF(:,92:5041); % trimmed for fair comparison
        timePoints = 92:5041;
    end
elseif strcmp(whichDepVar,'cumsumvel')
    depVar = cumsum(vel);
    if pooledPrediction|stackedPrediction
        depVar = depVar(92:5041);
        dFF = dFF(:,92:5041); % trimmed for fair comparison
        timePoints = 92:5041;
    end
elseif strcmp(whichDepVar,'displBout')
    depVar = displacementPerBout;
    dFF = dFF(:,92:5041); % trimmed for fair comparison
    timePoints = 92:5041;
    
elseif strcmp(whichDepVar,'boutType')
    boutType = calculateZebrafishBoutType(92,5041,vel);
    boutTypeCat = categorical(boutType);
    boutTypeOneHot = onehotencode(boutTypeCat,2);
    depVar = boutTypeOneHot;
    dFF = dFF(:,92:5041); % trimmed for fair comparison
    timePoints = 92:5041;
% elseif strcmp(whichDepVar,'displBoutAppropriate') %is the fish appropriately compensating? ideally, 0 displacement
%     depVar = calculateZebrafishDisplacementAppropriateness(displacementPerBout);
%     dFF = dFF(:,92:5041); % trimmed for fair comparison
%     timePoints = 92:5041;
elseif strcmp(whichDepVar,'time')
    depVar = timePoints.';
    if pooledPrediction|stackedPrediction
        depVar = depVar(92:5041);
        dFF = dFF(:,92:5041); % trimmed for fair comparison
        timePoints = 92:5041;
    end
elseif strcmp(whichDepVar,'bhvStdPerBout')
    bhvTrimmed = bhv(92:5041);
    E=reshape(bhvTrimmed,150,[]);
    nBouts = size(E,1);
    F=repmat(std(E,1),[1 150]);
    stdPerBout=reshape(F,[],1);
%     disp(size(stdPerBout))
    depVar = stdPerBout;

    dFF = dFF(:,92:5041); % trimmed for fair comparison
    timePoints = 92:5041;
elseif strcmp(whichDepVar,'displBoutPercent')
    depVar = 100.*displacementPerBout./sum(vel,"all");
    dFF = dFF(:,92:5041); % trimmed for fair comparison
    timePoints = 92:5041;
elseif strcmp(whichDepVar,'dvdt')
    depVar = [0; diff(vel)];
    if pooledPrediction|stackedPrediction
        depVar = depVar(92:5041);
        dFF = dFF(:,92:5041); % trimmed for fair comparison
        timePoints = 92:5041;
    end

elseif strcmp(whichDepVar,'dbdt')
    depVar = [0; diff(bhv)];
    if pooledPrediction|stackedPrediction
        depVar = depVar(92:5041);
        dFF = dFF(:,92:5041); % trimmed for fair comparison
        timePoints = 92:5041;
    end
    
%     dFF = dFF(:,1:end-1); % trimmed for fair comparison
%     timePoints = timePoints(1:end-1);
else
    depVar = bhv;
    if pooledPrediction|stackedPrediction
        depVar = depVar(92:5041);
        dFF = dFF(:,92:5041); % trimmed for fair comparison
        timePoints = 92:5041;
    end
    disp('default dependent variable is bhv')
end

lateStartShift = 842; %after first 5 complete bouts

if strcmp(whichPoints,'later') & (strcmp(whichDepVar,'displBout')|strcmp(whichDepVar,'boutType')|strcmp(whichDepVar,'bhvStdPerBout')|pooledPrediction|stackedPrediction|strcmp(whichDepVar,'displBoutPercent'))
    disp('here')
    depVar = depVar((lateStartShift-91):end,:);
    dFF = dFF(:,(lateStartShift-91):end);
    timePoints = timePoints((lateStartShift-91):end);
elseif strcmp(whichPoints,'later') 
    depVar = depVar(lateStartShift:end,:); 
    dFF = dFF(:,lateStartShift:end);
    timePoints = timePoints(lateStartShift:end);
end

end
%% map hd5
function varargout = mapFishHindbrainData(dataToMap,whichFilteredIdxLB)
% data to map could be variance explained per neuron, etc. - 1 dimensional
% whichFilteredIdxLB based on Ljung-Box test - use lbPerNeuronOriginal
% only these indices remain - must rematch these to original traces

volumeID = h5read("FishBehaviorData\20190923_ramp_bc_6hz\cells0_clean.hdf5", ...
    '/volume_id');
volumeMean = h5read("FishBehaviorData\20190923_ramp_bc_6hz\volume0.hdf5", ...
    '/volume_mean');

brainMap = permute(volumeMean,[3 2 1]);
brainMapMax = squeeze(max(brainMap,[],1));
% %make all non-zero cell values approximately same - no bright spots, just a
% %background
% brainMapNonzero    = brainMapMax;
% brainMapNonzero(brainMapNonzero == 0) = [];
% brainMapThreshold  = median(brainMapNonzero,"all");
% brainMapCapped     = brainMapMax;
% brainMapCapped(brainMapCapped>brainMapThreshold) = brainMapThreshold;
% brainMapCappedNotInf    = brainMapCapped;
% brainMapCappedNotInfIdx = isfinite(log(brainMapCapped));
% brainMapCappedNotInf(brainMapCappedNotInfIdx==0) = [];
% brainMapLogShift        = min(log(brainMapCappedNotInf),[],"all");
% brainMapLogShifted      = log(brainMapCapped)-brainMapLogShift;

% https://www.mathworks.com/matlabcentral/answers/476715-superimposing-two-imagesc-graphs-over-each-other
cla();

f1 = figure();
% %plot first data 
% ax1 = axes; 
% 
% im1 = imagesc(ax1,edge(im2gray(brainMapLogShifted)));
% 
% % im1 = imagesc(ax1,brainMapLogShifted);
% % clim([0 2*max(brainMapLogShifted,[],"all")])
% % colorbar
% % pbaspect([1 1 1])
% 
% im1.AlphaData = 0.3; % change this value to change the background image transparency 
% axis square; 
% hold all; 
% colormap(ax1,'summer') 

%plot second data 
ax2 = axes; 


% same trace used multiple times across locations on the map, goes from 0 
% to 14089 (originally python syntax, so added 1 to go from 1 to 14090
% instead)
%for whichSlice = 1:size(volumeID,1)
cellIDs     = squeeze(max(volumeID))+1;%squeeze(volumeID(whichSlice,:,:));%squeeze(max(volumeID))+1; %good enough? the largest possible cell id for that xy index
remappedIDs = 1:size(whichFilteredIdxLB,2);
remappedIDs(whichFilteredIdxLB(1,:)==0) = [];

% check - should be no Infs
if sum(isnan(dataToMap),'all') > 0
    warning('NaNs remain - check dataToMap');
end
% check - should be no Infs
if sum(isinf(dataToMap),'all') > 0
    warning('Infs remain - check dataToMap');
end

dataThresholdMax = max(dataToMap,[],'all');
dataThresholdMin = min(dataToMap,[],'all'); % FIXME: assumed

% brainMapOverlay      = -1000.*ones(size(brainMapMax));
brainMapOverlay      = nan(size(brainMapMax));
% brainMapOverlay      = zeros(size(brainMapMax));


% for i = 1:size(cellIDs,1)
    % for j = 1:size(cellIDs,2) % larger here
for i = 1:size(remappedIDs,2)
    [row,col]            =  find(cellIDs == remappedIDs(i));
    if ~isempty(row)
        % hold on
        for j = 1:size(row,1)
            brainMapOverlay(row(j),col(j)) = dataToMap(i); 
        end
        
    end
        
end



hold on
im2                      = imagesc(ax2,brainMapOverlay.');
% im2.AlphaData            = 0.5; % change this value to change the foreground image transparency 
% axis square; 
% axis([0 560 0 512])
% colormap(ax2,'winter') 
% clim([dataThresholdMin dataThresholdMax]);
axis xy
colorbar
hold off

brainMapOverlay(isnan(brainMapOverlay)) = 0;
% disp(sum(sum(brainMapOverlay(:,1:280))))
% disp(sum(sum(brainMapOverlay(:,281:560))))
disp(sum(sum(brainMapOverlay(:,1:270))))
disp(sum(sum(brainMapOverlay(:,271:512))))

disp((sum(sum(brainMapOverlay(:,1:270)))/sum(sum(~isnan(brainMapOverlay(:,1:270))))))
disp((sum(sum(brainMapOverlay(:,271:512)))/sum(sum(~isnan(brainMapOverlay(:,271:512))))))

% disp(sum(sum(dataToMap))) % not redundant even though the mapping is

% %link axes 
% linkaxes([ax1,ax2]) 
% %%Hide the top axes 
% ax2.Visible = 'off'; 
% ax2.XTick = []; 
% ax2.YTick = []; 
% %set the axes and colorbar position 
% set([ax1,ax2],'Position',[.17 .11 .685 .815]); 
% cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]); 
% cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]); 

% pause();
% cla(ax2);

%end
% varargout{1} = f1;

% filename = sprintf('varExplPerNeuronMap_En_%d.svg',);
% print('-painters','-dsvg',filename)

end
%% map hd5
function varargout = mapFishHindbrainConnectivity(dataToMap,whichFilteredIdxLB,Phi,...
    F,c_values,doTime,varargin)
% data to map could be variance explained per neuron, etc. - 1 dimensional
% whichFilteredIdxLB based on Ljung-Box test - use lbPerNeuronOriginal
% only these indices remain - must rematch these to original traces

% set(gca,"NextPlot","replacechildren")
v = VideoWriter("fishEn_test.avi");
v.FrameRate = 6; %6Hz
f1 = figure();
open(v)




doBlocks    = 0;
doCentroids = 0;
doActivity  = 0;
trialID = 'defaulttrial';

if nargin == 10
    traceVsBlock = varargin{1};
    blocksize1   = varargin{2};
    blocksize2   = varargin{3};
    doBlocks     = 1;
    trialID         = varargin{4};
elseif nargin == 9
    traceVsBlock = varargin{1};
    blocksize1   = varargin{2};
    blocksize2   = varargin{3};
    doBlocks     = 1;
elseif nargin == 8
    traceVsCentroid = varargin{1};
    doCentroids     = 1;
    trialID         = varargin{2};
elseif nargin == 7
    traceVsCentroid = varargin{1};
    doCentroids     = 1;
else
    error('Not set up yet - either input centroid locations or blocks and sizes')
end
    

volumeID = h5read("FishBehaviorData\20190923_ramp_bc_6hz\cells0_clean.hdf5", ...
    '/volume_id');
volumeMean = h5read("FishBehaviorData\20190923_ramp_bc_6hz\volume0.hdf5", ...
    '/volume_mean');



% for whichF = [5,17,27,29,33,36,48]%1:size(F,2)%19
    % unitnormF           = F{1,whichF}.*1./max(abs(F{1,whichF}),[],'all');
    % connectivityDueToFi = Phi*unitnormF*Phi.';
if doTime
    forLimit = size(c_values,1);
else
    forLimit = size(c_values,2);
    % doActivity = 1;
    % warning('Simulating activity due to each operator')
end

for snapshot = [8 10 18 25 2 5 16 24]%[8 10 18 25]%1:forLimit %1040:10:1180%1:forLimit
    % title("")
    if doTime
        timePoint      = snapshot;
        scaleEachDOByC = F;
        for whichDO = 1:size(F,2)
            scaleEachDOByC{whichDO} = c_values(timePoint,whichDO).*F{whichDO};
        end
        scaleEachDOByCArray = cat(3,scaleEachDOByC{:});
        dynsTogether        = sum(scaleEachDOByCArray,3);
        [U,L]               = eig(dynsTogether);
        maxEvalScaleFactor  = max(abs(L(:)));
        scaledF             = real(U*L*(diag(1./maxEvalScaleFactor)*inv(U)));
        connectivityDueToFi = Phi*scaledF*Phi.';
    else
        whichF              = snapshot;
        [U,L]               = eig(F{1,whichF});
        maxEvalScaleFactor  = max(abs(L(:)));
        unitnormF           = real(U*L*(diag(1./maxEvalScaleFactor)*inv(U)));
        connectivityDueToFi = Phi*unitnormF*Phi.';
    end
    
    % figure();histogram(connectivityDueToFi(:))
    % figure();imagesc(connectivityDueToFi);colorbar

    if ~isnan(std(connectivityDueToFi(:)))

        % impulse             = ones(size(F{1,whichF},1),1);
        % simxvals            = zeros(size(F{1,whichF},1),300);
        % simxvals(:,1)       = unitnormF * impulse;
        % for simTP = 2:300
        %     simxvals(:,simTP)  = unitnormF * simxvals(:,simTP-1);
        % end
        % activityDueToFi     = Phi*mean(simxvals,2,'omitmissing');
    
        % dataToMap = activityDueToFi;
    
        brainMap = permute(volumeMean,[3 2 1]);
        brainMapMax = squeeze(max(brainMap,[],1));
        
        % https://www.mathworks.com/matlabcentral/answers/476715-superimposing-two-imagesc-graphs-over-each-other
        
        %plot second data 
        ax2 = axes; 
        
        % same trace used multiple times across locations on the map, goes from 0 
        % to 14089 (originally python syntax, so added 1 to go from 1 to 14090
        % instead)
        %for whichSlice = 1:size(volumeID,1)
        cellIDs     = squeeze(max(volumeID))+1;%squeeze(volumeID(whichSlice,:,:));%squeeze(max(volumeID))+1; %good enough? the largest possible cell id for that xy index
        remappedIDs = 1:size(whichFilteredIdxLB,2);
        remappedIDs(whichFilteredIdxLB(1,:)==0) = [];
        
        % check - should be no Infs
        if sum(isnan(dataToMap),'all') > 0
            warning('NaNs remain - check dataToMap');
        end
        % check - should be no Infs
        if sum(isinf(dataToMap),'all') > 0
            warning('Infs remain - check dataToMap');
        end
        
        % dataThresholdMax = max(dataToMap,[],'all');
        % dataThresholdMin = min(dataToMap,[],'all'); % FIXME: assumed
        
        brainMapOverlay      = nan(size(brainMapMax));
        
        
        % for i = 1:size(cellIDs,1)
            % for j = 1:size(cellIDs,2) % larger here
        for i = 1:size(remappedIDs,2)
            [row,col]            =  find(cellIDs == remappedIDs(i));
            if ~isempty(row)
                % hold on
                for j = 1:size(row,1)
                    brainMapOverlay(row(j),col(j)) = dataToMap(i); 
                end
                
            end
                
        end
        
        hold on
        im2                      = imagesc(ax2,brainMapOverlay.');
        im2.AlphaData            = 0.2; % change this value to change the foreground image transparency 
        % axis([0 size(brainMapOverlay,2) 0 size(brainMapOverlay,1)])
        axis square; 
        % colormap(ax2,'winter') 
        % clim([dataThresholdMin dataThresholdMax]);
        % colorbar
        
        
        % connectivity
    
    
        if doBlocks
    
            title(sprintf("Operator %d",snapshot))
            centroidOfEachBlockX = int32(blocksize1/2):blocksize1:size(cellIDs,1); %approximate
            centroidOfEachBlockY = int32(blocksize2/2):blocksize2:size(cellIDs,2);
        
            connectionsBlockToBlockSum = zeros(int32(size(cellIDs,1)/blocksize1),...
                int32(size(cellIDs,2)/blocksize2),...
                int32(size(cellIDs,1)/blocksize1),...
                int32(size(cellIDs,2)/blocksize2));
        
            for mm = 1:int32(size(cellIDs,1)/blocksize1)
                for nn = 1:int32(size(cellIDs,2)/blocksize2)
                    whichTraces = find(traceVsBlock(:,mm,nn)==1);
                    if ~isempty(whichTraces)
                        sumConnectionsFromThisBlock = ...
                            sum(connectivityDueToFi(whichTraces,:),1); % sum for the target traces (down the columns) - should leave you with 1 row
        
                        for oo = 1:int32(size(cellIDs,1)/blocksize1)
                            for pp = 1:int32(size(cellIDs,2)/blocksize2)
                                % whichTraces2 = find(traceVsBlock(:,oo,pp)==1);
                                connectionsBlockToBlockSum(mm,nn,oo,pp) = ...
                                    sum(sumConnectionsFromThisBlock(1,traceVsBlock(:,oo,pp)==1),"all"); % one number
                            end
                        end
                        
                    end
                end
            end
        
            rescaleFactor = max(abs(connectionsBlockToBlockSum),[],'all');
            thresholdToPlot = 0.3; %10*std(connectionsBlockToBlockSum(:));
            disp(thresholdToPlot)
            disp(sum(find(abs(connectionsBlockToBlockSum)>thresholdToPlot)))
            connectionsBlockToBlockSum(abs(connectionsBlockToBlockSum)/rescaleFactor < thresholdToPlot) = 0;
            connectionsBlockToBlockSum = connectionsBlockToBlockSum/rescaleFactor;
        
            for mm = 1:int32(size(cellIDs,1)/blocksize1)
                for nn = 1:int32(size(cellIDs,2)/blocksize2)
                    for oo = 1:int32(size(cellIDs,1)/blocksize1)
                        for pp = 1:int32(size(cellIDs,2)/blocksize2)
                            
        
                            p1 = [centroidOfEachBlockX(mm) centroidOfEachBlockY(nn)];   % First Point
                            p2 = [centroidOfEachBlockX(oo) centroidOfEachBlockY(pp)];       % Second Point
                            dp = p2-p1;                         % Difference
        
                            if connectionsBlockToBlockSum(mm,nn,oo,pp) ~= 0
                                if connectionsBlockToBlockSum(mm,nn,oo,pp) > 0
                                    arrowColor = 'r'; % hot
                                else
                                    arrowColor = 'b'; % cold
                                end
                                
                                hold on
                                if sum(dp,'all') ~= 0
                                    fprintf('\n%d %d %d %d %.2f', mm,nn,oo,pp,abs(connectionsBlockToBlockSum(mm,nn,oo,pp))/rescaleFactor)
                                    quiver(p1(1),p1(2),dp(1),dp(2),...
                                        'LineWidth',(abs(connectionsBlockToBlockSum(mm,nn,oo,pp))),...
                                        'ShowArrowHead','on',...
                                        'Color',arrowColor)
                                end
                            end
                        end
                    end
                    
                    
                end
            end
    
        elseif doCentroids

            if doActivity
                tic        
                impulseFOutput = ones(size(Phi,2),500);
                for genPoint = 2:500
                    impulseFOutput(:,genPoint) = unitnormF*impulseFOutput(:,genPoint-1);
                end
                yImpulse         = Phi * impulseFOutput;
                rescaledActivity = yImpulse/max(abs(yImpulse),[],'all');
                rescaledActPos   = rescaledActivity;
                rescaledActNeg   = rescaledActivity;
                rescaledActPos(rescaledActivity<0) = NaN;
                rescaledActNeg(rescaledActivity>0) = NaN;
                rescaledActNeg   = abs(rescaledActNeg);
                
                % activityDueToFi = zeros(size(cellIDs,1),size(cellIDs,2));
                for genPoint = 1:500
                    title(sprintf("Operator %d, simulated time point %d",snapshot,genPoint))
                    hold on
                    im2                      = imagesc(ax2,brainMapOverlay.');
                    im2.AlphaData            = 0.2; % change this value to change the foreground image transparency 
                    % axis([0 size(brainMapOverlay,2) 0 size(brainMapOverlay,1)])
                    axis square; 
                    scatter(traceVsCentroid(:,2),traceVsCentroid(:,3),1000*rescaledActPos(:,genPoint),'red') % hot
                    scatter(traceVsCentroid(:,2),traceVsCentroid(:,3),1000*rescaledActNeg(:,genPoint),'blue') % cold
                    toc
                    axis([0 size(brainMapOverlay,1) 0 size(brainMapOverlay,2)])
                    axis square
                    hold off
                    % varargout{1} = f1;
                    frame = getframe(gcf);
                    writeVideo(v,frame)
                    cla
                end
                clf reset
            else
                rescaleFactor = max(abs(connectivityDueToFi),[],'all');
                thresholdToPlot = 0.1; %10*std(connectionsBlockToBlockSum(:));
                connectivityDueToFi(abs(connectivityDueToFi)/rescaleFactor < thresholdToPlot) = 0;
                connectivityDueToFi = connectivityDueToFi/rescaleFactor;


                tic
                if doTime
                    title(sprintf("Frame %d of %d",snapshot,size(c_values,1)))
                else
                    title(sprintf("Operator %d of %d",snapshot,size(c_values,2)))
                end
                for ww = 1:size(remappedIDs,2)
                    % traceIDSource   = traceVsCentroid(ww,1);
                    centroidXSource = traceVsCentroid(ww,2);
                    centroidYSource = traceVsCentroid(ww,3);
                    
                    if ~isnan(centroidXSource) & ~isnan(centroidYSource)
        
                        for yy = 1:size(remappedIDs,2)
                            
                            % traceIDTarget   = traceVsCentroid(yy,1);
                            centroidXTarget = traceVsCentroid(yy,2);
                            centroidYTarget = traceVsCentroid(yy,3);
        
                            if ~isnan(centroidXTarget) & ~isnan(centroidYTarget)
                                % thresholdToPlot = 1e-2;
                                p1 = [centroidXSource centroidYSource];   % First Point
                                p2 = [centroidXTarget centroidYTarget];   % Second Point
                                dp = p2-p1;                               % Difference
                    
                                if connectivityDueToFi(ww,yy) ~= 0
                                    if connectivityDueToFi(ww,yy) > 0
                                        arrowColor = 'r'; % hot
                                    else
                                        arrowColor = 'b'; % cold
                                    end
                        
                                    hold on
                                    if sum(dp,'all') ~= 0
                                        fprintf('\n%d %d %d %d %.2f', centroidXSource,centroidYSource,centroidXTarget,centroidYTarget,abs(connectivityDueToFi(ww,yy)))
                                        quiver(p1(1),p1(2),dp(1),dp(2),...
                                            'LineWidth',abs(connectivityDueToFi(ww,yy)),...
                                            'ShowArrowHead','on',...
                                            'Color',arrowColor)
                                    end
                                end
                            
                            end
                        end
                    end
                end
                toc
            end
        end
        
        
        fprintf('\n');
    
        if ~doActivity
            axis([0 size(brainMapOverlay,1) 0 size(brainMapOverlay,2)])
            axis square
            hold off
            % varargout{1} = f1;

            % trialID = input('TrialID:','s');
            % trialID = "bhv_250428_lambdap1";
            if doTime
                filename = sprintf("%s_Frame%d.svg",trialID,snapshot);
            else
                filename = sprintf("%s_DO%d.svg",trialID,snapshot);
            end
            print('-painters','-dsvg',filename)

            frame = getframe(gcf);
            writeVideo(v,frame)
            % pause();
        end
    end
    % cla;
    clf reset
end
close(v)

end

%%
function traceVsBlock = lookupTraceVsBlock(whichFilteredIdxLB,bs1,bs2)

volumeID = h5read("FishBehaviorData\20190923_ramp_bc_6hz\cells0_clean.hdf5", ...
    '/volume_id');

cellIDs     = squeeze(max(volumeID))+1;%squeeze(volumeID(whichSlice,:,:));%squeeze(max(volumeID))+1; %good enough? the largest possible cell id for that xy index
remappedIDs = 1:size(whichFilteredIdxLB,2);
remappedIDs(whichFilteredIdxLB(1,:)==0) = [];

% store which block(s) each remappedID belongs to
blocksize1 = bs1; % divisible into 512 and 560 
blocksize2 = bs2; % divisible into 512 and 560 
traceVsBlock = zeros(size(remappedIDs,2),...
    int32(size(cellIDs,1)/blocksize1),...
    int32(size(cellIDs,2)/blocksize2));

for kk = 1:size(remappedIDs,2)
    disp(kk)
    for mm = 1:int32(size(cellIDs,1)/blocksize1)
        for nn = 1:int32(size(cellIDs,2)/blocksize2)
            relevantCellIdx1 = blocksize1*(mm-1)+1:blocksize1*(mm-1)+blocksize1;
            relevantCellIdx2 = blocksize2*(nn-1)+1:blocksize2*(nn-1)+blocksize2;
            valuesInBlock = unique(cellIDs(relevantCellIdx1,relevantCellIdx2));
            % disp(valuesInBlock)
            % disp(remappedIDs(kk))
            if ~isempty(find(valuesInBlock==remappedIDs(kk),1))
                traceVsBlock(kk,mm,nn) = 1;
            end
        end
    end
end

end
%%
function traceVsCentroid = lookupTraceVsCentroid(whichFilteredIdxLB)

volumeID = h5read("FishBehaviorData\20190923_ramp_bc_6hz\cells0_clean.hdf5", ...
    '/volume_id');

cellIDs     = squeeze(max(volumeID))+1;%squeeze(volumeID(whichSlice,:,:));%squeeze(max(volumeID))+1; %good enough? the largest possible cell id for that xy index
remappedIDs = 1:size(whichFilteredIdxLB,2);
remappedIDs(whichFilteredIdxLB(1,:)==0) = [];

% store which block(s) each remappedID belongs to
traceVsCentroid = zeros(size(remappedIDs,2),3);

for kk = 1:size(remappedIDs,2)
    disp(kk)
    traceVsCentroid(kk,1) = remappedIDs(kk);
    [idxX,idxY]           = find(cellIDs==remappedIDs(kk));
    idxCentrX             = round(mean(idxX));
    idxCentrY             = round(mean(idxY));
    traceVsCentroid(kk,2) = idxCentrX;
    traceVsCentroid(kk,3) = idxCentrY;
end

end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%