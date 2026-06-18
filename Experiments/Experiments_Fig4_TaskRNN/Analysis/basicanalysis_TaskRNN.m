%%
clear
%%
addpath(genpath('.'))
%%
bdlds_01 = load('saveTaskRNN_Psiregbdlds_260423iters500bdlds_s04.mat');
% bdlds_01 = load('./TaskRNNData/saveTaskRNN_afterBO_260422opt4params_s01_extravarsdeleted.mat');
% bdlds_02 = load('./TaskRNNData/saveTaskRNN_afterBO_260422opt4params_s02.mat');
% bdlds_03 = load('./TaskRNNData/saveTaskRNN_afterBO_260422opt4params_s03.mat');
%%
bdldsx_01 = load("saveTaskRNN_Psireg_260429_s01.mat");
% bdldsx_01 = load("saveTaskRNN_Psireg_260429_s00_vemultiopt_v2.mat")
% bdldsx_01 = load("saveTaskRNN_Psireg_260427veBhvOnly_s05.mat");
% bdldsx_01 = load("saveTaskRNN_Psireg_260427_s00_vemultiopt.mat")
% bdldsx_01 = load('./TaskRNNData/saveTaskRNN_afterBO_260423opt5paramsbdldsx_s01.mat');
% bdldsx_02 = load('./TaskRNNData/saveTaskRNN_afterBO_260423opt5paramsbdldsx_s02.mat');
% bdldsx_03 = load('./TaskRNNData/saveTaskRNN_afterBO_260423opt5paramsbdldsx_s03.mat');
%%
% meanBhvR2bdlds = mean([bdlds_01.varExplBhv,bdlds_02.varExplBhv,bdlds_03.varExplBhv])
% %%
% meanBhvR2bdldsx = mean([bdldsx_01.varExplBhv,bdldsx_02.varExplBhv,bdldsx_03.varExplBhv])

%% check outputs - b-dLDS
figure()
tiledlayout('vertical')
nexttile()
plot(bdlds_01.dFF{1}.');ylabel('activations')
box off
nexttile()
plot(bdlds_01.behaviorData{1}.');ylabel('bhv')
box off
nexttile()
plot(bdlds_01.A_cell{1}.');ylabel('x')
box off
nexttile()
plot(bdlds_01.B_cell{1}.');ylabel('c')
box off
nexttile()
plot((bdlds_01.Psi*bdlds_01.B_cell{1}).');ylabel('Psi*c');subtitle(sprintf('Behavior R2 %0.2f',bdlds_01.varExplBhv))
box off
%% check outputs - b-dLDS-x
figure()
tiledlayout('vertical')
nexttile()
plot(bdldsx_01.dFF{1}.');ylabel('activations')
box off
nexttile()
plot(bdldsx_01.behaviorData{1}.');ylabel('bhv')
box off
nexttile()
plot(bdldsx_01.A_cell{1}.');ylabel('x')
box off
nexttile()
plot(bdldsx_01.B_cell{1}.');ylabel('c')
box off
nexttile()
plot((bdldsx_01.Psi*bdldsx_01.A_cell{1}).');ylabel('Psi*x');subtitle(sprintf('Behavior R2 %0.2f',bdldsx_01.varExplBhv))
box off
%% use behavior R^2 to quantify success (A_cell is data, so neural reconstruction R^2 not super relevant here)
% reconstructedBehavior = Psi*A_cell{1}; % for b-dLDS-x
% allbhv = behaviorData{1};
% rval = corrcoef(allbhv(:), reconstructedBehavior(:));
% varExplBhv = (rval(1,2))^2;
%%
whichModel = bdlds_01;

varExplBhvPerX = zeros(size(whichModel.A_cell{1},1),size(whichModel.behaviorData{1},1));
for i = 1:size(whichModel.A_cell{1},1)
    for j = 1:size(whichModel.behaviorData{1},1)
        rBhvX             = corrcoef(whichModel.A_cell{1}(i,:),whichModel.behaviorData{1}(j,:));
        varExplBhvPerX(i,j) = rBhvX(1,2).^2;
    end
end

varExplBhvPerC = zeros(size(whichModel.B_cell{1},1),size(whichModel.behaviorData{1},1));
for i = 1:size(whichModel.B_cell{1},1)
    for j = 1:size(whichModel.behaviorData{1},1)
        rBhvX             = corrcoef(whichModel.B_cell{1}(i,:),whichModel.behaviorData{1}(j,:));
        varExplBhvPerC(i,j) = rBhvX(1,2).^2;
    end
end
%%
figure();imagesc(bdlds_01.Psi);colorbar;title('Psi, b-dLDS');ylabel('bhv');xlabel('c');box off; grid off
%%
figure();imagesc(bdldsx_01.Psi);colorbar;title('Psi, b-dLDS-x');ylabel('bhv');xlabel('x');box off; grid off
%%
[maxVEx,idxmaxVEx] = max(varExplBhvPerX(:))
[xidx, bhvidx] = ind2sub(size(varExplBhvPerX), idxmaxVEx)
[maxVEc,idxmaxVEc] = max(varExplBhvPerC(:))
[cidx, bhvidx] = ind2sub(size(varExplBhvPerC), idxmaxVEc)
%%
[maxPsi, idxMaxPsi] = max(whichModel.Psi(:))
[bhvidx, mappedidx] = ind2sub(size(whichModel.Psi(:)), idxMaxPsi)
%%
medVEx = median(varExplBhvPerX(:))
stdVEx = std(varExplBhvPerX(:))
medVEc = median(varExplBhvPerC(:))
stdVEc = std(varExplBhvPerC(:))
%%
medPsi = median(whichModel.Psi(:))
stdPsi = std(whichModel.Psi(:))
%%
figure();imagesc(bdlds_01.Psi);colorbar;title('Psi, b-dLDS, lambda=200');xlabel('c');ylabel('behaviors');box off
%%
figure();imagesc(bdldsx_01.Psi);colorbar;title('Psi, b-dLDS-x, lambda=200');xlabel('x');ylabel('behaviors');box off
%% L1/L2 norm (sum abs)/sum(squares)
whichModel = bdlds_01;
L1L2 = sum(abs(whichModel.Psi(:)))/sum(whichModel.Psi.^2,"all")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%