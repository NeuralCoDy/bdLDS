addpath(genpath('.'))
%%
clear; close all; clc
%%
load("saveFMT_1bhvgraded_forCLDS_GOOD_251113.mat")
load("C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/clds/notebooks/cldsoutput.mat")
%% construct 8x8x3000(x50) for learned dLDS model (lin. comb. f and c)
fcMatDLDS = zeros(8,8,3000,50); % one time trace for each entry in combined Fc matrix, for each sample
for i = 1:length(B_cell) % samples N_ex
    for j = 1:length(F) % learned operators
        thisC  = B_cell{i}(j,:); % 1 by 3000 time points
        thisF  = F{j};
        % thisFC = repmat(thisC,[8 8 1]);
        for k = 1:3000 % time points
            fcMatDLDS(:,:,k,i) = fcMatDLDS(:,:,k,i) + thisF*thisC(k);
        end
    end
end
%%
cldsPermute  = permute(clds,[2 3 1]);
%%
tiledlayout(3,3)
nexttile()
imagesc(trueFC(:,:,750,1))
colorbar
clim([-1.5 1.5])
grid off
colormap redbluecmap
box off
title('Frame 750')
ylabel('Ground truth')
nexttile()
imagesc(trueFC(:,:,1250,1))
colorbar
clim([-1.5 1.5])
grid off
colormap redbluecmap
box off
title('Frame 1250')
nexttile()
imagesc(trueFC(:,:,2250,1))
colorbar
clim([-1.5 1.5])
grid off
colormap redbluecmap
box off
title('Frame 2250')


nexttile()
imagesc(fcMatDLDS(:,:,750,1))
colorbar
clim([-1.5 1.5])
grid off
colormap redbluecmap
box off
ylabel('b-dLDS')
nexttile()
imagesc(fcMatDLDS(:,:,1250,1))
colorbar
clim([-1.5 1.5])
grid off
colormap redbluecmap
box off
nexttile()
imagesc(fcMatDLDS(:,:,2250,1))
colorbar
clim([-1.5 1.5])
grid off
colormap redbluecmap
box off

nexttile()
imagesc(cldsPermute(:,:,750))
colorbar
clim([-1.5 1.5])
grid off
colormap redbluecmap
box off
ylabel('CLDS')
nexttile()
imagesc(cldsPermute(:,:,1250))
colorbar
clim([-1.5 1.5])
grid off
colormap redbluecmap
box off
nexttile()
imagesc(cldsPermute(:,:,2250))
colorbar
clim([-1.5 1.5])
grid off
colormap redbluecmap
box off
%% CLDS neural reconstruction: x_t+1 = A_t * x_t
x_sample = dFF{1}; % x = y (D=I)
nDim = size(x_sample,1);
nTp  = size(x_sample,2);
x_t_hat = zeros(nDim,nTp-1);
for tp = 2:nTp
    x_t_hat(:,tp) = squeeze(clds(tp,:,:))*x_sample(:,tp-1);
end

figure(), cla;
tSel = 1; % which sample
nLines = 2;
subplot(nLines,4,[1,4]), plot(dFF{tSel}(1,:).')       % ASC added 7/25
box off; xlabel('Time (frames)'); ylabel('Signal (GT)')
set(gca,'XLim',[1,size(dFF{tSel},2)]);
subplot(nLines,4,[5,8]), plot(x_t_hat(1,:).')    % ASC added 7/25
set(gca,'XLim',[1,size(dFF{tSel},2)]);
box off; xlabel('Time (frames)'); ylabel('CLDS Reconstruction')


%%
% cm  = [1 0 0; 1 1 1; 0 0 1];                     % Basic Colormap
% cmi = interp1([-1.5; 0; 1.5], cm, (-1.5:1.5));         % interpolated Colormap
%% plot top left corner (16 rows) - 3 heat maps (true, dlds, clds)
% topleftTrueFC = reshape(trueFC(1:4,1:4,:,1),16,3000); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
% topleftDLDSFC = reshape(fcMatDLDS(1:4,1:4,:,1),16,3000); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
% topleftCLDSFC = reshape(cldsPermute(1:4,1:4,:),16,3000); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
% fig1 = figure();
% tiledlayout('vertical')
% nexttile()
% imagesc(topleftTrueFC)
% colorbar
% title('True At=sum(Fm*cmt), rows 1:4, cols 1:4 (top left corner)')
% yticks([1 5 9 13])
% yticklabels({'A11', 'A12', 'A13', 'A14'})
% % yticklabels({'r1c1', 'r2c1', 'r3c1', 'r4c1',...
% %     'r1c2', 'r2c2', 'r3c2', 'r4c2',...
% %     'r1c3', 'r2c3', 'r3c3', 'r4c3',...
% %     'r1c4', 'r2c4', 'r3c4', 'r4c4'})
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
% 
% nexttile()
% imagesc(topleftDLDSFC)
% colorbar
% title('dLDS At=sum(Fm*cmt) (learned), rows 1:4, cols 1:4 (top left corner)')
% yticks([1 5 9 13])
% yticklabels({'A11', 'A12', 'A13', 'A14'})
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
% 
% nexttile()
% imagesc(topleftCLDSFC)
% colorbar
% title('CLDS Dynamics matrix At, rows 1:4, cols 1:4 (top left corner)')
% yticks([1 5 9 13])
% yticklabels({'A11', 'A12', 'A13', 'A14'})
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
% 
% %% plot top right corner (16 rows) - 3 heat maps (true, dlds, clds)
% toprightTrueFC = reshape(trueFC(1:4,5:8,:,1),16,3000); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
% toprightDLDSFC = reshape(fcMatDLDS(1:4,5:8,:,1),16,3000); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
% toprightCLDSFC = reshape(cldsPermute(1:4,5:8,:),16,3000); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
% fig2 = figure();
% tiledlayout('vertical')
% nexttile()
% imagesc(toprightTrueFC)
% colorbar
% title('True At=sum(Fm*cmt), rows 1:4, cols 5:8 (top right corner)')
% yticks([1 5 9 13])
% yticklabels({'A15', 'A16', 'A17', 'A18'})
% % yticklabels({'r1c1', 'r2c1', 'r3c1', 'r4c1',...
% %     'r1c2', 'r2c2', 'r3c2', 'r4c2',...
% %     'r1c3', 'r2c3', 'r3c3', 'r4c3',...
% %     'r1c4', 'r2c4', 'r3c4', 'r4c4'})
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
% 
% nexttile()
% imagesc(toprightDLDSFC)
% colorbar
% title('dLDS At=sum(Fm*cmt) (learned), rows 1:4, cols 5:8 (top right corner)')
% yticks([1 5 9 13])
% yticklabels({'A15', 'A16', 'A17', 'A18'})
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
% 
% nexttile()
% imagesc(toprightCLDSFC)
% colorbar
% title('CLDS Dynamics matrix At, rows 1:4, cols 5:8 (top right corner)')
% yticks([1 5 9 13])
% yticklabels({'A15', 'A16', 'A17', 'A18'})
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
% 
% %% plot bottom left corner (16 rows) - 3 heat maps (true, dlds, clds)
% bottomleftTrueFC = reshape(trueFC(5:8,1:4,:,1),16,3000); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
% bottomleftDLDSFC = reshape(fcMatDLDS(5:8,1:4,:,1),16,3000); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
% bottomleftCLDSFC = reshape(cldsPermute(5:8,1:4,:),16,3000); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
% fig3 = figure();
% tiledlayout('vertical')
% nexttile()
% imagesc(bottomleftTrueFC)
% colorbar
% title('True At=sum(Fm*cmt), rows 5:8, cols 1:4 (bottom left corner)')
% yticks([1 5 9 13])
% yticklabels({'A51', 'A52', 'A53', 'A54'})
% % yticklabels({'r1c1', 'r2c1', 'r3c1', 'r4c1',...
% %     'r1c2', 'r2c2', 'r3c2', 'r4c2',...
% %     'r1c3', 'r2c3', 'r3c3', 'r4c3',...
% %     'r1c4', 'r2c4', 'r3c4', 'r4c4'})
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
% 
% nexttile()
% imagesc(bottomleftDLDSFC)
% colorbar
% title('dLDS At=sum(Fm*cmt) (learned), rows 5:8, cols 1:4 (bottom left corner)')
% yticks([1 5 9 13])
% yticklabels({'A51', 'A52', 'A53', 'A54'})
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
% 
% nexttile()
% imagesc(bottomleftCLDSFC)
% colorbar
% title('CLDS Dynamics matrix At, rows 5:8, cols 1:4 (bottom left corner)')
% yticks([1 5 9 13])
% yticklabels({'A51', 'A52', 'A53', 'A54'})
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
% 
% %% plot bottom right corner (16 rows) - 3 heat maps (true, dlds, clds)
% bottomrightTrueFC = reshape(trueFC(5:8,5:8,:,1),16,3000); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
% bottomrightDLDSFC = reshape(fcMatDLDS(5:8,5:8,:,1),16,3000); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
% bottomrightCLDSFC = reshape(cldsPermute(5:8,5:8,:),16,3000); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
% fig4 = figure();
% tiledlayout('vertical')
% nexttile()
% imagesc(bottomrightTrueFC)
% colorbar
% title('True At=sum(Fm*cmt), rows 5:8, cols 5:8 (bottom right corner)')
% yticks([1 5 9 13])
% yticklabels({'A55', 'A56', 'A57', 'A58'})
% % yticklabels({'r1c1', 'r2c1', 'r3c1', 'r4c1',...
% %     'r1c2', 'r2c2', 'r3c2', 'r4c2',...
% %     'r1c3', 'r2c3', 'r3c3', 'r4c3',...
% %     'r1c4', 'r2c4', 'r3c4', 'r4c4'})
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
% 
% nexttile()
% imagesc(bottomrightDLDSFC)
% colorbar
% title('dLDS At=sum(Fm*cmt) (learned), rows 5:8, cols 5:8 (bottom right corner)')
% yticks([1 5 9 13])
% yticklabels({'A55', 'A56', 'A57', 'A58'})
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
% 
% nexttile()
% imagesc(bottomrightCLDSFC)
% colorbar
% title('CLDS Dynamics matrix At, rows 5:8, cols 5:8 (bottom right corner)')
% yticks([1 5 9 13])
% yticklabels({'A55', 'A56', 'A57', 'A58'})
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
