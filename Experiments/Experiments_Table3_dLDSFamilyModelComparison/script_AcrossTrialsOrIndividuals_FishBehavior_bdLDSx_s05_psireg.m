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
save('saveFish_En_260416_bdLDSx_s05_psireg.mat','-v7.3');
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