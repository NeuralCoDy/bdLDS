%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

% parpool(16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code
addpath(genpath('.'))
%%
rng('shuffle')
currentState1 = rng;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate data
% use 'sizeD' latent states, 'sizeD' "recorded channels", 'sizeT' timepoints

sizeD            = 8;
sizeT            = 3000; % T/2
simulatedDmatrix = eye(sizeD,sizeD);
contOpt          = false;
nF               = [3,3];



[latentStatesX, Fgt, groundTruthStates] = generateContMultiSubNetwork(...
               sizeT,[sizeD/2,sizeD/2],nF,50,contOpt,false,false,true);
%%

nBhv = 10;

% fourierPsi = fourierBasis(sum(nF,"all"));
% rowPsi = fourierPsi(rand);
% rowPsi = rowPsi(1,sum(nF,"all")+1:end);
% simulatedPsi = repmat(rowPsi,nBhv,1);

% %%
% A = rand(nBhv,sum(nF,'all'));
% [nrows, ncols] = size(A);
% Q = zeros(nrows, ncols);
% R = zeros(ncols,ncols);
% for j = 1:ncols
%     v = A(:,j);
%     for i = 1:j-1
%         R(i,j) = Q(:,i).'*A(:,j);
%         v = v-R(i,j)*Q(:,i);
%     end
%     R(j,j) = norm(v);
%     Q(:,j) = v/R(j,j);
% end
% simulatedPsi = Q;

% nCoefs = sum(nF,'all');
% fourierPsi = fourierBasis(nBhv);
% % randomFrequencies = randsample(nCoefs,nCoefs);
% % randomOffset = rand(nCoefs,1)*pi/2;
% simulatedPsi = zeros(nBhv,nCoefs);
% 
% for i = 1:sum(nF,'all')
%     colPsi = (fourierPsi(i)).';
%     colPsi = colPsi(nBhv+1:end,1);
%     simulatedPsi(:,i) = colPsi;
% end
% 
% simulatedPsi = simulatedPsi + rand(nBhv,nCoefs);


A = randn(nBhv,sum(nF,'all')); % sum(nF,'all')

% A(:,2:end) = 0; 

% %try 2 coeffs tied to behavior
% A(:,3:end) = 0; % first 2 coefficients are related to behavior
% A(round(nBhv/2)+1:end,1) = 0; %dyn coeff 1 generates behaviors 1-5
% A(1:round(nBhv/2),2) = 0; % dyn coeff 2 generates behaviors 6-10

% try overlapping
% A(:,3:end) = 0; % first 2 coefficients are related to behavior
% A(round(nBhv/2)+2:end,1) = 0; %dyn coeff 1 generates behaviors 1-7
% A(1:round(nBhv/2)-2,2) = 0; % dyn coeff 2 generates behaviors 3-10
% 

simulatedPsi = A * 1./vecnorm(A);
simulatedPsi(isnan(simulatedPsi))=0;



% simulatedPsi = eye(nBhv);

sum(simulatedPsi.^2,1)
simulatedPsi.'*simulatedPsi
corr(simulatedPsi)
disp('check for correlation between columns - if too high, might impede reconstruction of behavior')
justOneGTState = groundTruthStates;
% for cellID = 1:size(groundTruthStates,1)
%     justOneGTState{cellID,1}(2:end,:) = 0;
% end
behaviorData = cellfun(@(x) simulatedPsi*x,justOneGTState,'UniformOutput',false);
% disp('FIXME: can behaviorData be multiple cells too?')
%%

%try more noise (last arg: usually 0)
dFF = statesToObservations(simulatedDmatrix, latentStatesX, 0);
% figure();cla;plot(dFF{1}.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bayesopt
% addpath(genpath('~/my_documents/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code'))

% rng default
% lambda_behavior = optimizableVariable('lambda_behavior',[0.001 5],'Type','real'); % [0.001 1] --> 0.5536
% % lambda_b = optimizableVariable('lambda_b',[0.01 2],'Type','real'); % [0.001 1] --> 0.5536
% % lambda_historyb = optimizableVariable('lambda_historyb',[0.01 2],'Type','real'); % [0.1 10], real --> 9.9994
% % fun = @(x)runBehaviorDLDSSimLambdaHistoryBLambdaBehavior(dFF,behaviorData,...
% %     latentStatesX,groundTruthStates,simulatedPsi,x.lambda_behavior,x.lambda_historyb);
% % fun = @(x)runBehaviorDLDSSim_LambdaBLambdaHistoryB_BestDynCoeffs(dFF,behaviorData,...
%     % latentStatesX,groundTruthStates,simulatedPsi,x.lambda_b,x.lambda_historyb);
% fun = @(x)runBehaviorDLDSSimLambdaBehavior(dFF,behaviorData,...
%     latentStatesX,groundTruthStates,simulatedPsi,x.lambda_behavior);
% results = bayesopt(fun,lambda_behavior,'Verbose',0,...
%     'AcquisitionFunctionName','expected-improvement-plus');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Behavior and Psi - see after dLDS
% 
% % Can coefficients be learned from behavior?
% 
% whichSample = 1;
% 
% for eachCoeff = 1:size(groundTruthStates{whichSample},1)
%     predictorB = behaviorData{whichSample}.';
%     outputOneC = groundTruthStates{whichSample}(eachCoeff,:).';
%     [weights,fitinfo] = lasso(predictorB,outputOneC);
%     [~,idxMSE] = min(fitinfo.MSE);
%     figure()
%     hold on
%     plot(outputOneC)
%     plot(weights(:,idxMSE).'*predictorB.')
%     legend('actual c','predicted c')
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Make ground truth
% 
% Fgt2 = cell(numel(Fgt{1})+numel(Fgt{2})+1,1);
% for ll = 1:numel(Fgt2);   Fgt2{ll} = zeros(sizeD,sizeD);              end
% for ll = 1:numel(Fgt{1}); Fgt2{ll}(1:sizeD/2,1:sizeD/2) = Fgt{1}{ll}; end
% for ll = 1:numel(Fgt{2})
%     Fgt2{ll+numel(Fgt{1})}(sizeD/2+1:sizeD,sizeD/2+1:sizeD) = Fgt{2}{ll}; 
% end
% Fgt2{end}                        = 0.001*randn(sizeD,sizeD);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Set parameters
inf_opts.nF              = 15; % # dynamics operators
inf_opts.N               = size(latentStatesX{1},1); % # latent states
inf_opts.M               = size(dFF{1},1) ; % original dimension (# channels)
% disp('If reset dFF below, make sure to change inf_opts.M');
inf_opts.lambda_val      = 0.0001;         % ASC added 7/25 % VARY - tradeoff: approx. SNR - see ratios of traces % 0.0001 (don't shrink this further - gets to be essentially 0 for solver_L1RLS)
inf_opts.lambda_history  = 0.0001;%0.01;         % ASC added 7/25
inf_opts.lambda_b        = input('lambda_b:');%0.25; %0.5;%0.39;           % ASC added 7/25 %0.025 % 0.5 %0.6
inf_opts.lambda_historyb = input('lambda_historyb:');%1.16;           % ASC added 7/25 %0.7 %0.45 %0.35 % 1.99 way too strong - all coefficients end up same
inf_opts.tol             = 1e-8;           % 1e-3
inf_opts.max_iter2       = 20;             % 500 %20
inf_opts.max_iters       = input('max iters:');
inf_opts.F_update        = true; % default = true;
inf_opts.D_update        = false; % NoObs case - default = true;
inf_opts.N_ex            = input('N_ex (e.g., 50):');%50;   
inf_opts.T_s             = input('T_s (e.g., 200):');%200; 
if inf_opts.T_s > size(dFF{1},2)
    warning('Sample duration (time points T_s) too large for data')
    fprintf('\nDefaulting to sample size == size of data %d\n',size(dFF{1},2))
    inf_opts.T_s = size(dFF{1},2);
end
inf_opts.step_d          = 1;  % ASC added 7/25
inf_opts.step_f          = 10;  % ASC added 7/25 % 30
inf_opts.step_decay      = input('step_decay (e.g., 0.9999):'); %0.998 %0.9999 helped get closer - 4/6 operators in corners 
inf_opts.plot_option     = 10;  % ASC added 7/25
inf_opts.lambda_f_decay  = input('lambda_f_decay (e.g., 0.9999):'); %0.996 %0.996
inf_opts.lambda_f        = input('lambda_f (e.g., 0.2):'); %0.2 almost there but two operators nearly identical %0.4 too high - only 1-2 kinds of operators emerge
inf_opts.solver_type     = input('Solver (fista or tfocs)?','s');%'tfocs'; % fista
inf_opts.special         = 'noobs';
inf_opts.deltaDynamics   = false; %default: false %use x_t-x_{t-1}
inf_opts.AcrossIndividuals     = 0; %input('Across individuals? 1 for yes, 0 for no:');

inf_opts.behaviordLDS    = input('behavior_dLDS?:');
if inf_opts.behaviordLDS == 1
    inf_opts.step_psi        = input('step_psi:');%1; 
    inf_opts.lambda_behavior = input('lambda_behavior:'); %1.6; %0.5521; % 0.5536 for opt for psi
    inf_opts.verysparsebhv   = 1;
    inf_opts.lambda2         = input('lambda2:');
    inf_opts.psinorm         = input('Norm on Psi - norm or frob?:','s');
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

% [Phi, F] = bpdndf_dynamics_learning(dFF, [], simulatedDmatrix, inf_opts);              % Run the dictionary learning algorithm
% % Phi = simulatedDmatrix; F = Fgt2;
% [A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFF, Phi, F, ...
%                                        @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
%%
% addpath(genpath('~/my_documents/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/'))
%%
rng('shuffle')
currentState2 = rng;


if inf_opts.behaviordLDS
    [Phi, F, Psi] = bpdndf_dynamics_learning_behavior(dFF, [], [], [], behaviorData, inf_opts);              % Run the dictionary learning algorithm
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
%%

if inf_opts.behaviordLDS
    for whichSample = 1:size(B_cell,1)
        plotMultiNetworkOutput(whichSample, dFF, A_cell, B_cell, groundTruthStates, F, Psi, simulatedPsi,behaviorData)
        plotBehaviorOnly(whichSample, dFF, A_cell, B_cell, groundTruthStates, F, Psi, simulatedPsi,behaviorData)
        pause()
        % close all
    end
else
    for whichSample = 1:size(B_cell,1)
        plotMultiNetworkOutput(whichSample, dFF, A_cell, B_cell, groundTruthStates, F)
        pause()
        % close all
    end
end

%% Can coefficients be learned from behavior?

% whichSample = 1;
% 
% for eachCoeff = 1:size(B_cell{whichSample},1)
%     predictorB = behaviorData{whichSample}.';
%     outputOneC = B_cell{whichSample}(eachCoeff,:).';
%     [weights,fitinfo] = lasso(predictorB,outputOneC);
%     [~,idxMSE] = min(fitinfo.MSE);
%     figure()
%     hold on
%     plot(outputOneC)
%     plot(weights(:,idxMSE).'*predictorB.')
%     legend('actual c','predicted c')
% end

%%
varExplX = varExplLatentsOnly(A_cell,latentStatesX);
varExplX_summary = mean(varExplX);
%%
% varExplC = varExplLatentsOnly(B_cell,groundTruthStates);

% Since there are different numbers of ground truth and inferred coeffs,
% each combination gets evaluated and returned for each sample. Then, we
% evaluate the maximum variance explained for each ground truth coeff among
% the inferred coeffs. Then we take the mean of these across all of the
% states and samples (in a cell array).
% disp("omitmissing is more recent - use omitnan for r2022a")
% varExplC_summary = mean(cellfun(@(x) mean(max(x,[],1),"all","omitnan"),varExplC), "all"); % should be close to 1
% 

%% compare Psi to ground truth Psi

varExplPsi = varExplLatentsOnly({simulatedPsi.'},{Psi.'});
disp("omitmissing is more recent - use omitnan for r2022a")
varExplPsi_summary = mean(cellfun(@(x) mean(max(x,[],2),"all","omitnan"),varExplPsi), "all"); % should be close to 1

%%
% figure()
% set(gcf,"Color","white")
% imagesc(simulatedPsi)
% colorbar
% title('Ground-truth Psi')
% 
% %%
% figure()
% set(gcf,"Color","white")
% imagesc(Psi)
% colorbar
% title('Learned Psi')
% 
%%
for whichSample = 1:size(behaviorData,1)
    figure()
    set(gcf,"Color","white")
    subplot(2,4,1:4)
    if exist('Psi','var')
        cmapMaxDistinct = distinguishable_colors(size(Psi,1));
    else
        cmapMaxDistinct = distinguishable_colors(size(simulatedPsi,1));
    end
    % cmapMaxDistinct = distinguishable_colors(size(Psi,1));
    hold on
    tVals = 1:size(behaviorData{whichSample},2);
    for ii=1:size(behaviorData{whichSample},1)
        Y = behaviorData{whichSample};
        plot(tVals,Y(ii,:),'color',cmapMaxDistinct(ii,:),'DisplayName',num2str(ii));
        hold on;
    end
    legend
    ylabel('True behavior')
    hold off

    subplot(2,4,5:8)
    hold on
    for ii=1:size(Psi,1)
        Y = (Psi*B_cell{whichSample});
        plot(tVals,Y(ii,:),'color',cmapMaxDistinct(ii,:),'DisplayName',num2str(ii));
        hold on;
    end
    legend
    ylabel('Reconstructed behavior')
    xlabel('Time')
    hold off
    pause();
end
%%
% % figure();
% varExplPerBehaviorTrace = zeros(size(behaviorData,1),size(behaviorData{1},1));
% for i = 1:size(behaviorData,1) % samples
%     for j = 1:size(behaviorData{1},1) % behavior dimensions
%         behaviorIJ = behaviorData{i}(j,:); % one behavior trace by time, from one sample
%         reconstrBehaviorI = Psi*B_cell{i}; % all behavior traces by time, from one sample
%         reconstrBehaviorIJ = reconstrBehaviorI(j,:);
%         % % find spots where 1 frame increases, next frame decreases, or vice versa (unnatural jump)
%         % idxIncrease = diff(reconstrBehaviorIJ)>1*std(reconstrBehaviorIJ);
%         % idxDecrease = diff(reconstrBehaviorIJ)<-1*std(reconstrBehaviorIJ);
%         % 
%         % behaviorIJ(find(idxIncrease)) = behaviorIJ(find(idxIncrease)+1);
%         % reconstrBehaviorIJ(find(idxIncrease)) = reconstrBehaviorIJ(find(idxIncrease)+1);
%         % behaviorIJ(find(idxDecrease)) = behaviorIJ(find(idxDecrease)+1);
%         % reconstrBehaviorIJ(find(idxDecrease)) = reconstrBehaviorIJ(find(idxDecrease)+1);
%         % 
%         % hold on
%         % plot(behaviorIJ.')
%         % plot(reconstrBehaviorIJ.')
%         % legend('True behavior','Reconstructed behavior')
%         % hold off
%         % pause(1)
%         % cla
%         r = corrcoef(behaviorIJ,reconstrBehaviorIJ);
%         varExplPerBehaviorTrace(i,j) = (r(1,2))^2;
%     end
% end
% 
% varExplBhv = mean(varExplPerBehaviorTrace,'all');
% %%
% figure()
% set(gcf,"Color","white")
% histogram(varExplPerBehaviorTrace,50)
% title('Var Expl per behavior trace, all samples')
% 
% % %%
% % r = corrcoef(behaviorData{1}.',(Psi*B_cell{1}).');
% % varExplBhv = (r(1,2))^2;
% %%
% figure()
% tiledlayout("flow")
% for i = 1:size(F,2)
%     nexttile
%     imagesc(F{i})
% end
% 
% %%
% % h5disp('cells0_clean.hdf5')
% cell_x = h5read('cells0_clean.hdf5','/cell_x');
% cell_y = h5read('cells0_clean.hdf5','/cell_y');
% cell_z = h5read('cells0_clean.hdf5','/cell_z');
% %neuronLocs = 
% hold on
% for i = 1:size(F,2) 
%     figure()
%     connections = Phi*F{i}*Phi.';
%     impulse = ones(size(Phi,2),1);
%     x_justthisF = zeros(size(Phi,2),300);
%     for l = 1:300
%         x_justthisF(:,l) = F{i}*impulse;
%     end
%     y_justThisF = Phi*x_justthisF;
%     ymean = mean(y_justThisF,2); % each trace
%     ymean01 = (ymean-min(ymean))/(max(ymean)-min(ymean));
%     for j = 1:size(connections,1)
%         hold on
%         scatter3(vj(1), vj(2), vj(3), ymean01(j))
%         for k = 1:size(connections,2)
%             hold on
%             vj=[j j j]; % fix these with neuron locations
%             vk=[k k k];
%             v=[vj;vk];
%             plot3(v(:,1),v(:,2),v(:,3),'r')
%             hold off
%         end
%     end
% 
% end
% hold off
%% recreate BayesOpt landscape if didn't save
% xA = results.XTrace.lambda_behavior;
% yA = results.XTrace.lambda_historyb;
% zA = results.ObjectiveMinimumTrace;
% N = 250;
% xvec = linspace(min(xA), max(xA), N);
% yvec = linspace(min(yA), max(yA), N);
% [X, Y] = ndgrid(xvec, yvec);
% F = scatteredInterpolant(xA, yA, zA);
% Z = F(X, Y);
% figure()
% surf(X, Y, Z, 'edgecolor', 'none');
% xlabel('lambda b')
% ylabel('lambda history b')
% figure()
% scatter3(xA, yA, zA, [], zA);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%