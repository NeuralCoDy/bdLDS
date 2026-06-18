load('saveFMT_1bhvgraded_forCLDS_TMLR40_rerun_GOOD.mat')
%%
if inf_opts.behaviordLDS
    for whichSample = 1:size(B_cell,1)
        % plotMultiNetworkOutput(whichSample, dFF, A_cell, B_cell, groundTruthStates, F, Psi, simulatedPsi)
        plotMultiNetworkOutput(whichSample, dFF, A_cell, B_cell, groundTruthStates, F, Psi, simulatedPsi,behaviorData)
        plotBehaviorOnly(whichSample, dFF, A_cell, B_cell, groundTruthStates, F, Psi, simulatedPsi,behaviorData,1)
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
varExplX = varExplLatentsOnly(A_cell(11:50),latentStatesX(11:50));
varExplX_summary = mean(varExplX);
%%
varExplC = varExplLatentsOnly(B_cell(11:50),groundTruthStates(11:50));

% Since there are different numbers of ground truth and inferred coeffs,
% each combination gets evaluated and returned for each sample. Then, we
% evaluate the maximum variance explained for each ground truth coeff among
% the inferred coeffs. Then we take the mean of these across all of the
% states and samples (in a cell array).
% disp("omitmissing is more recent - use omitnan for r2022a")
varExplC_summary = mean(cellfun(@(x) mean(max(x,[],1),"all","omitnan"),varExplC), "all"); % should be close to 1


%% compare Psi to ground truth Psi

% varExplPsi = varExplLatentsOnly({Psi.'},{simulatedPsi.'});
% disp("omitmissing is more recent - use omitnan for r2022a")
% varExplPsi_summary = mean(cellfun(@(x) mean(max(x,[],2),"all","omitnan"),varExplPsi), "all"); % should be close to 1
%%
Fgt2 = cell(numel(Fgt{1})+numel(Fgt{2}),1);
for ll = 1:numel(Fgt2);   Fgt2{ll} = zeros(sizeD,sizeD);              end
for ll = 1:numel(Fgt{1}); Fgt2{ll}(1:sizeD/2,1:sizeD/2) = Fgt{1}{ll}; end
for ll = 1:numel(Fgt{2})
    Fgt2{ll+numel(Fgt{1})}(sizeD/2+1:sizeD,sizeD/2+1:sizeD) = Fgt{2}{ll}; 
end


varExplF = varExplLatentsOnlyFComparison(Fgt2,F.');
figure();imagesc(varExplF);xlabel('Inferred dyn. ops.');ylabel('Ground truth dyn. ops.');title('R^2')
varExplF_summary = mean(max(varExplF,[],2)); %max of each row (ground truth state reconstruction)

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
% for whichSample = 1:size(behaviorData,1)
%     figure()
%     set(gcf,"Color","white")
%     subplot(2,4,1:4)
%     if exist('Psi','var')
%         cmapMaxDistinct = distinguishable_colors(size(Psi,1));
%     else
%         cmapMaxDistinct = distinguishable_colors(size(simulatedPsi,1));
%     end
%     % cmapMaxDistinct = distinguishable_colors(size(Psi,1));
%     hold on
%     tVals = 1:size(behaviorData{whichSample},2);
%     for ii=1:size(behaviorData{whichSample},1)
%         Y = behaviorData{whichSample};
%         plot(tVals,Y(ii,:),'color',cmapMaxDistinct(ii,:),'DisplayName',num2str(ii));
%         hold on;
%     end
%     legend
%     ylabel('True behavior')
%     hold off
% 
%     subplot(2,4,5:8)
%     hold on
%     for ii=1:size(Psi,1)
%         Y = (Psi*B_cell{whichSample});
%         plot(tVals,Y(ii,:),'color',cmapMaxDistinct(ii,:),'DisplayName',num2str(ii));
%         hold on;
%     end
%     legend
%     ylabel('Reconstructed behavior')
%     xlabel('Time')
%     hold off
%     pause();
% end
%%
% % figure();
% varExplPerBehaviorTrace = zeros(size(behaviorData,1),size(behaviorData{1},1));
varExplPerBehaviorTrace = zeros(50,size(behaviorData{1},1));
for i = 11:50%size(behaviorData,1) % samples
    for j = 1:size(behaviorData{1},1) % behavior dimensions
        behaviorIJ = behaviorData{i}(j,:); % one behavior trace by time, from one sample
        reconstrBehaviorI = Psi*B_cell{i}; % all behavior traces by time, from one sample
        reconstrBehaviorIJ = reconstrBehaviorI(j,:);
        % % find spots where 1 frame increases, next frame decreases, or vice versa (unnatural jump)
        % idxIncrease = diff(reconstrBehaviorIJ)>1*std(reconstrBehaviorIJ);
        % idxDecrease = diff(reconstrBehaviorIJ)<-1*std(reconstrBehaviorIJ);
        % 
        % behaviorIJ(find(idxIncrease)) = behaviorIJ(find(idxIncrease)+1);
        % reconstrBehaviorIJ(find(idxIncrease)) = reconstrBehaviorIJ(find(idxIncrease)+1);
        % behaviorIJ(find(idxDecrease)) = behaviorIJ(find(idxDecrease)+1);
        % reconstrBehaviorIJ(find(idxDecrease)) = reconstrBehaviorIJ(find(idxDecrease)+1);
        % 
        % hold on
        % plot(behaviorIJ.')
        % plot(reconstrBehaviorIJ.')
        % legend('True behavior','Reconstructed behavior')
        % hold off
        % pause(1)
        % cla
        r = corrcoef(behaviorIJ,reconstrBehaviorIJ);
        varExplPerBehaviorTrace(i,j) = (r(1,2))^2;
    end
end

varExplBhv = mean(varExplPerBehaviorTrace(11:50,:),'all','omitnan');