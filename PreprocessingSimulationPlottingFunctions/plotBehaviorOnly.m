function plotBehaviorOnly(tSel, dFF, A_cell, B_cell, groundTruthStates, F,varargin)
if nargin > 9
    nLines = 7;%9;
    Psi = varargin{1};
    simulatedPsi = varargin{2};
    behaviorData = varargin{3};
    ifPlotBothOnOnePlot = varargin{4};
elseif nargin > 8
    nLines = 7;%9;
    Psi = varargin{1};
    simulatedPsi = varargin{2};
    behaviorData = varargin{3};
    ifPlotBothOnOnePlot = 0;
elseif nargin > 6
    nLines = 7;
    Psi = varargin{1};
    simulatedPsi = varargin{2};
    ifPlotBothOnOnePlot = 0;
else
    nLines = 5;
    ifPlotBothOnOnePlot = 0;
end


% tSel = 5; % selected trial
figure(), cla;
% subplot(nLines,4,[1,4]), plot(dFF{tSel}.')       % ASC added 7/25
% box off; xlabel('Time (frames)'); ylabel('Signal (GT)')
% set(gca,'XLim',[1,size(dFF{tSel},2)]);
% subplot(nLines,4,[5,8]), plot(A_cell{tSel}.')    % ASC added 7/25
% set(gca,'XLim',[1,size(dFF{tSel},2)]);
% box off; xlabel('Time (frames)'); ylabel('Reconstruction')
% 
% subplot(nLines,4,[9,12]), imagesc(abs(B_cell{tSel}))   % ASC added 7/25
% % subplot(nLines,4,[1,4]), imagesc(abs(B_cell{tSel}))   % ASC added 7/25
% 
% set(gca,'XLim',[1,size(dFF{tSel},2)]);
% box off; xlabel('Time (frames)'); ylabel('Coeff Index')
% 
% subplot(5,4,[9,12]), plot(B_cell{tSel}.')   % ASC added 7/25
% set(gca,'XLim',[1,size(dFF{tSel},2)]);
% box off; xlabel('Time (frames)'); ylabel('Dynamics coefficient amplitude')
% legend('C1','C2','C3','C4')
% 
% subplot(nLines,4,[13,16]), imagesc(groundTruthStates{tSel})   % ASC added 7/25
% % subplot(nLines,4,[5,8]), imagesc(groundTruthStates{tSel})   % ASC added 7/25
% 
% box off; xlabel('Time (frames)'); ylabel('True state'); ylim('padded')
% %set(gca,'XLim',[1,size(dFF{tSel},2)],'YLim',[-0.1,1.1]);
% % legend('C1','C2','C3','C4')
% 
% allFs = F{1};
% nF1   = size(F{1},1);
% for ll = 2:numel(F)
%     allFs = cat(2,allFs,[zeros(nF1,5),F{ll}]);
% end
% subplot(nLines,4,[17,20]), imagesc(allFs); 
% % subplot(nLines,4,[9,12]), imagesc(allFs); 
% box off; ylabel('Learned DOs');
% axis image; axis off; colormap gray


% [~, maxind] = max(abs(B_cell{1}),[],2);

if nargin > 6
    % % subplot(nLines,4,[21,24]), imagesc(Psi); box off; ylabel('Learned Psi');
    % % subplot(nLines,4,[25,28]), imagesc(simulatedPsi); box off; ylabel('True Psi');
    % 
    % disp('bhvall 250910 Psi')
    % 
    % 
    % hold on
    % % subplot(nLines,4,[21,24]); 
    % subplot(nLines,4,[13,16]); 
    % 
    % selectPsi    = [B_cell{tSel}(5,440)*Psi(:,5) ...
    %     B_cell{tSel}(12,440)*Psi(:,12) ...
    %     B_cell{tSel}(5,440)*Psi(:,5)+B_cell{tSel}(12,440)*Psi(:,12)...
    %     B_cell{tSel}(6,196)*Psi(:,6) ...
    %     B_cell{tSel}(10,196)*Psi(:,10) ...
    %     B_cell{tSel}(6,196)*Psi(:,6)+B_cell{tSel}(10,196)*Psi(:,10) ...
    %     B_cell{tSel}(4,3000)*Psi(:,4) ...
    %     B_cell{tSel}(1,500)*Psi(:,1) ...
    %     B_cell{tSel}(8,1081)*Psi(:,8) ...
    %     B_cell{tSel}(13,2241)*Psi(:,13)].';   
    % % selectPsi    = [B_cell{tSel}(5,maxind(5))*Psi(:,5) ...
    % %     B_cell{tSel}(12,maxind(12))*Psi(:,12) ...
    % %     B_cell{tSel}(5,maxind(5))*Psi(:,5)+B_cell{tSel}(12,maxind(12))*Psi(:,12)...
    % %     B_cell{tSel}(6,maxind(6))*Psi(:,6) ...
    % %     B_cell{tSel}(10,maxind(10))*Psi(:,10) ...
    % %     B_cell{tSel}(6,maxind(6))*Psi(:,6)+B_cell{tSel}(10,maxind(10))*Psi(:,10) ...
    % %     B_cell{tSel}(4,maxind(4))*Psi(:,4) ...
    % %     B_cell{tSel}(1,maxind(1))*Psi(:,1) ...
    % %     B_cell{tSel}(8,maxind(8))*Psi(:,8) ...
    % %     B_cell{tSel}(13,maxind(13))*Psi(:,13)].';    
    % % selectPsi    = [B_cell{tSel}(5,maxind(5))*Psi(:,5) ...
    % %     B_cell{tSel}(12,maxind(12))*Psi(:,12) ...
    % %     B_cell{tSel}(5,maxind(5))*Psi(:,5)+B_cell{tSel}(12,maxind(12))*Psi(:,12)...
    % %     B_cell{tSel}(6,maxind(6))*Psi(:,6) ...
    % %     B_cell{tSel}(10,maxind(10))*Psi(:,10) ...
    % %     B_cell{tSel}(6,maxind(6))*Psi(:,6)+B_cell{tSel}(10,maxind(10))*Psi(:,10)].';   
    % bar(1:10,selectPsi) %,'DisplayName','Learned Psi'
    % set(gca,'XTick',1:numel(1:10), 'XTickLabel',['5   '; '12  '; '5&12'; '6   '; '10  '; '6&10'; '4   '; '1   '; '8   '; '13  '])
    % % bar(1:6,selectPsi) %,'DisplayName','Learned Psi'
    % % set(gca,'XTick',1:numel(1:10), 'XTickLabel',['1'; '1'; '1'; '2'; '2'; '2'])
    % % xtickangle(90); 
    % ylabel('Learned Psi')
    % 
    % % selectPsi    = B_cell{tSel}(5,maxind(5))*Psi(:,5);           
    % % bar(selectPsi) %,'DisplayName','Learned Psi'
    % 
    % 
    % % bar(1:8,[Psi(:,5 12 6 10 4 1 8 13]); simulatedPsi(:,[1 1 2 2 3 4 5 6])],'grouped') %,'DisplayName','Learned Psi'
    % % bar(1:8,simulatedPsi(:,[1 1 2 2 3 4 5 6]),'grouped','DisplayName','True Psi')
    % box off;
    % hold off
    % 
    % hold on
    % % subplot(nLines,4,[25,28]); 
    % subplot(nLines,4,[17,20]); 
    % 
    % selectSimPsi = [Psi(:,1) Psi(:,1) Psi(:,1) Psi(:,2) Psi(:,2) Psi(:,2) Psi(:,3) Psi(:,4) Psi(:,5) Psi(:,6)].';
    % % selectSimPsi = [Psi(:,1) Psi(:,1) Psi(:,1) Psi(:,2) Psi(:,2) Psi(:,2)].';
    % bar(1:10,selectSimPsi) %,'DisplayName','Learned Psi'
    % set(gca,'XTick',1:numel(1:10), 'XTickLabel',['1'; '1'; '1'; '2'; '2'; '2'; '3'; '4'; '5'; '6'])
    % % bar(1:6,selectSimPsi) %,'DisplayName','Learned Psi'
    % % set(gca,'XTick',1:numel(1:10), 'XTickLabel',['1'; '1'; '1'; '2'; '2'; '2'])
    % % xtickangle(90); 
    % ylabel('True Psi')
    % 
    % % bar(1:8,[Psi(:,5 12 6 10 4 1 8 13]); simulatedPsi(:,[1 1 2 2 3 4 5 6])],'grouped') %,'DisplayName','Learned Psi'
    % % bar(1:8,simulatedPsi(:,[1 1 2 2 3 4 5 6]),'grouped','DisplayName','True Psi')
    % box off;
    % hold off

    if ifPlotBothOnOnePlot
        hold on
        % subplot(nLines,4,[21,24])
        plot(behaviorData{tSel}.','DisplayName','True behavior')       % ASC added 7/25
        % subplot(nLines,4,[29,32]), plot(behaviorData{tSel}.')       % ASC added 7/25
        box off; xlabel('Time (frames)'); ylabel('Behavior')
        behaviorReco = Psi*B_cell{tSel};
        plot(behaviorReco.','DisplayName','Reconstructed behavior')
        legend
        set(gca,'XLim',[1,size(dFF{tSel},2)]);
    else
        subplot(nLines,4,[21,24]), plot(behaviorData{tSel}.')       % ASC added 7/25
        box off; xlabel('Time (frames)'); ylabel('Behavior')
        set(gca,'XLim',[1,size(dFF{tSel},2)]);
        behaviorReco = Psi*B_cell{tSel};
        subplot(nLines,4,[25,28]), plot(behaviorReco.')       % ASC added 7/25
        box off; xlabel('Time (frames)'); ylabel('bhv reco')
        set(gca,'XLim',[1,size(dFF{tSel},2)]);

    end
end

fontsize(scale=0.5)

% subplot(5,4,17), imagesc(F{1}); axis image; axis off; colormap gray
% subplot(5,4,18), imagesc(F{2}); axis image; axis off; colormap gray
% subplot(5,4,19), imagesc(F{3}); axis image; axis off; colormap gray
% subplot(5,4,20), imagesc(F{4}); axis image; axis off; colormap gray

end