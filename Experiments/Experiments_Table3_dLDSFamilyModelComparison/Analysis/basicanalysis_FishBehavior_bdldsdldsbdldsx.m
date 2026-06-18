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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
plotSummaryCandReconstr_OneX(dFF,Phi,B_cell,A_cell,inf_opts,trialIdentifier,varExplPerNeuron,samplingRate,dataIdentifier);

%%
varExpl = overallVarExpl(dFF,A_cell,Phi,inf_opts);
%%
varExplPerNeuron = perNeuronVarExpl(dFF,A_cell,Phi,inf_opts);
%% if b-dLDS-x
reconstructedBehavior = Psi*A_cell{1};
allbhv = bhvcell{1};
rval = corrcoef(allbhv(:), reconstructedBehavior(:));
varExplBhv = (rval(1,2))^2;
%% if b-dLDS
% bhvcell = behaviorData;
reconstructedBehavior = Psi*B_cell{1};
allbhv = bhvcell{1};
rval = corrcoef(allbhv(:), reconstructedBehavior(:));
varExplBhv = (rval(1,2))^2;

%%
varExplBhvPerX = zeros(size(A_cell{1},1),size(behavior_data{1},1));
for i = 1:size(A_cell{1},1)
    for j = 1:size(behavior_data{1},1)
        rBhvX             = corrcoef(A_cell{1}(i,:),behavior_data{1}(j,:));
        varExplBhvPerX(i,j) = rBhvX(1,2).^2;
    end
end

varExplBhvPerC = zeros(size(B_cell{1},1),size(behavior_data{1},1));
for i = 1:size(B_cell{1},1)
    for j = 1:size(behavior_data{1},1)
        rBhvX             = corrcoef(B_cell{1}(i,:),behavior_data{1}(j,:));
        varExplBhvPerC(i,j) = rBhvX(1,2).^2;
    end
end
%%
[maxVEx,idxmaxVEx] = max(varExplBhvPerX(:))
[maxVEc,idxmaxVEc] = max(varExplBhvPerC(:))
%%
[maxPsi, idxMaxPsi] = max(Psi(:))
%%
medVEx = median(varExplBhvPerX(:))
stdVEx = std(varExplBhvPerX(:))
medVEc = median(varExplBhvPerC(:))
stdVEc = std(varExplBhvPerC(:))
%%
medPsi = median(Psi(:))
stdPsi = std(Psi(:))
%%
% ifDisplayOnSecondMonitor = input('Display on second monitor?(1 yes, 0 no):');
% p = get(0, "MonitorPositions");
% c_values = B_cell{1}.';
% Phit = Phi.';
% 
% for timePoint = 200 %[200 400 600 800 1000 1200 1400]
%     f = figure();
%     if ifDisplayOnSecondMonitor
%         f.Position = p(2,:);
%     else
%         f.Position = p(1,:);
%     end
%     f.WindowState = 'maximized';
% 
%     scaleEachDOByC = F;
%     for whichDO = 1:size(F,2)
%         scaleEachDOByC{whichDO} = c_values(timePoint,whichDO).*F{whichDO};
%     end
%     scaleEachDOByCArray = cat(3,scaleEachDOByC{:});
%     dynsTogether =  sum(scaleEachDOByCArray,3);
%     [U,L]   = eig(dynsTogether);
%     maxEvalScaleFactor = max(abs(L(:)));
%     scaledF = real(U*L*(diag(1./maxEvalScaleFactor)*inv(U)));
% 
%     hold on;
%     connectivityMatrix                          = Phi * scaledF * Phit;
%     connThresholded                             = connectivityMatrix;
%     upperQuantileThresholdStrength              = quantile(abs(connectivityMatrix),0.9999,"all");
%     connThr                                     = upperQuantileThresholdStrength; %0.18 * max(abs(connectivityMatrix(:))); % 1 stdev
%     connThresholded(abs(connThresholded) < connThr) = 0;
%     connThresholdedAbs = abs(connThresholded);
%     circularGraph(connThresholdedAbs);
%     % circularGraph(connThresholdedAbs,...
%     %     'Label',cellstr(channelNamesStrings));
% 
%     hold off
% end


%%
%  find(abs(Psi) > 0.2)
%  find(abs(Psi) < 0.05)
% %%
% for whichBhv = 1:4
%     if whichBhv == 1
%         find(abs(Psi(whichBhv,:)) > 0.5)
%         % find(abs(Psi(whichBhv,:)) < 0.25)
%     else
%         find(abs(Psi(whichBhv,:)) > 0.95)
%         % find(abs(Psi(whichBhv,:)) < 0.01)
%     end
% end
% %%
% for whichBhv = 1:4
%     [maxvals, maxidx] = maxk(abs(Psi(whichBhv,:)),4)
% end
% %%
% figure()
% for i = 1:25 %[8 10 18 25 2 5 16 24]%1:size(B_cell{1})
%     for j = 1:inf_opts.nBhv
%         cla
% 
%         nTimepoints = size(A_cell{1},2);
%         tVals = linspace(0,nTimepoints/samplingRate,nTimepoints);
% 
%         % tiledlayout(4,1)
%         tiledlayout('vertical')
% 
%         nexttile
%         plot(tVals,B_cell{1}(i,:))
%         box off
%         axis([0 max(tVals) -Inf Inf]);
%         rBhvCo = corrcoef(B_cell{1}(i,:),behavior_data{1}(j,:));
%         title(sprintf('Dyn coef %d, Psi coef %0.2f, R2 vs behavior %d %0.2f',i,Psi(j,i), j, rBhvCo(1,2)^2))
% 
%         nexttile
%         % plot(behavior_data{1}(2,:))
%         plot(tVals,behavior_data{1})
%         box off
%         axis([0 max(tVals) -Inf Inf]);
%         title('Behavior (smoothed), input to model')
% 
% 
%         nexttile
%         bhvReconstruction = Psi*B_cell{1};
%         plot(tVals,bhvReconstruction)
%         box off
%         axis([0 max(tVals) -Inf Inf]);
%         title('Behavior reconstruction')
%         legend(['motor (1)'; '2        '; '3        '; '4        '])
%         xlabel('Time (s)')
% 
% 
%         if doSmoothedBehavior
%             sigma           = smoothingwindow; 
%             gaussian_range  = -3*sigma:3*sigma; % setting up Gaussian window
%             gaussian_kernel = normpdf(gaussian_range,0,sigma); % setting up Gaussian kernel
%             gaussian_kernel = gaussian_kernel/sum(gaussian_kernel);
%             depVarRe        = conv(bhvReconstruction(:),gaussian_kernel(:),'same');
% 
%         end
% 
%         % depVarRe = depVarRe./(max(abs(depVarRe))); %1D behavior
% 
%         % nexttile
%         % plot(tVals,depVarRe)
%         % box off
%         % axis([0 max(tVals) -Inf Inf]);
%         % title('Behavior reconstruction (smoothed)')
% 
% 
%         rBhvRe = corrcoef(depVarRe,behavior_data{1});
%         % rBhvCo = corrcoef(B_cell{1}(i,:),behavior_data{1}(2,:));
% 
%         disp(i)
%         disp(rBhvRe(1,2)^2)
%         disp(rBhvCo(1,2)^2)
%         disp(Psi(i))
% 
%         pause()
% 
%     end
% end
% %%
% whichData = varExplPerNeuron;
% histVarExpl(whichData,inf_opts,trialIdentifier,dataIdentifier);
% %%
% mapFishHindbrainData(varExplPerNeuron,lbPerNeuronOriginal);
% %%
% mapFishHindbrainData(dFF{1},lbPerNeuronOriginal);
% %%
% % traceVsBlock = lookupTraceVsBlock(lbPerNeuronOriginal,64,70);
% % traceVsBlock = lookupTraceVsBlock(lbPerNeuronOriginal,1,1);
% traceVsCentroid = lookupTraceVsCentroid(lbPerNeuronOriginal);
% 
% c_values = B_cell{1}.';
% %%
% % mapFishHindbrainConnectivity(varExplPerNeuron,lbPerNeuronOriginal,Phi,F,c_values,0,traceVsBlock,64,70)
% % mapFishHindbrainConnectivity(varExplPerNeuron,lbPerNeuronOriginal,Phi,F,c_values,0,traceVsBlock,1,1)
% %% operators
% mapFishHindbrainConnectivity(varExplPerNeuron,lbPerNeuronOriginal,Phi,F,c_values,0,traceVsCentroid,trialIdentifier); %do operators
% %% time points
% mapFishHindbrainConnectivity(varExplPerNeuron,lbPerNeuronOriginal,Phi,F,c_values,1,traceVsCentroid,trialIdentifier); %do timepoints
% 
% %%
% if inf_opts.AcrossIndividuals
%     for ii = 1:size(Phi,1)
%         x_values = cell2mat(A_cell{ii});
%         dataReconstructed = (Phi{ii} * x_values).';    
% 
%         flatDFF = reshape(dFF{ii}, numel(dFF{ii}),[]);
%         flatReconstructed = reshape(dataReconstructed,numel(dataReconstructed),[]);
%         rval = corrcoef(dFF{ii}, dataReconstructed.');
%         varExpl = (rval(1,2))^2;
%         disp(varExpl);
% 
%     end
% else
%     for ii = 1:size(A_cell,1) %in the event of multiple trials
%         x_values = A_cell{ii};
%         dataReconstructed = (Phi * x_values).';    
% 
%         maxvalDFF = max(dFF{ii},[],'all');
%         stdzdOn = false;
%         if stdzdOn
%             dataReconstructedRescaled = maxvalDFF .*(Phi * x_values);
%             dFFRescaled = maxvalDFF .* dFF{ii};
%         else
%             dFFRescaled = dFF{ii};
%             dataReconstructedRescaled = (Phi * x_values);
%         end
% 
% 
%         flatDFF = reshape(dFFRescaled, numel(dFFRescaled),[]);
%         flatReconstructed = reshape(dataReconstructedRescaled,numel(dataReconstructedRescaled),[]);
%         rval = corrcoef(flatDFF, flatReconstructed);
%         varExpl = (rval(1,2))^2;
%         disp(varExpl);
%         % pause();
%     end
% end
% %%
% [~,whichTracesMin] = mink(varExplPerNeuron(:,ii),1000);
% [~,whichTracesMax] = maxk(varExplPerNeuron(:,ii),1000);
% 
% %%
% figure()
% tiledlayout(4,1)
% nexttile
% imagesc(dFF{1}(whichTracesMin,:))
% title('orig worst')
% colorbar
% nexttile
% imagesc(dataReconstructedRescaled(whichTracesMin,:))
% colorbar
% title('recon worst')
% 
% nexttile
% imagesc(dFF{1}(whichTracesMax,:))
% colorbar
% title('orig best')
% nexttile
% imagesc(dataReconstructedRescaled(whichTracesMax,:))
% colorbar
% title('recon best')
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

for snapshot = [18 1 16 23 20 10 12 6 7 3 22 2 9 8 24 15]%[8 10 18 25 2 5 16 24]%[8 10 18 25]%1:forLimit %1040:10:1180%1:forLimit
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