%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parpool(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code
addpath(genpath('.'))
% whichModelFile = input('Which Emmanuel fish model?','s');
% load(whichModelFile)
%%
disp('Change path to modeEstimates for your setup')
ifForCIS = input('CIS?(1 or 0):');
if ifForCIS
    addpath(genpath('~/my_documents/CoDybase/CoDybase-MATLAB/CoDybase-MATLAB/stats/'))
else
    addpath(genpath('../../../CoDybase/CoDybase-MATLAB/CoDybase-MATLAB/'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% load('saveFish_EmNew1_251103.mat')
% load('saveFish_EmNew1_251103_Bcell.mat')
% load('C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/ZebrafishData/FishDataEmmanuel123/Fish1/Fish1_centroids_aligned_MF59_masks1st.mat')
% coords = my_data;
if ifForCIS
    % load('FishDataEmmanuel123/Fish1/Fish1cfos.mat')
    % load('FishDataEmmanuel123/Fish1/save_Emmanuel_NewFish1_dLDS.mat')
    load('FishDataEmmanuel123/Fish1/save_EmNewFish1_dLDS_250105.mat', 'B_cell','F','Phi','filtered','varExpl','varExplPerNeuron')
    load('FishDataEmmanuel123/Fish1/save_EmNewFish1_250109_allgenesidxonly.mat', 'filtered')
    disp('Reminder: used dFF, not filtered.dFF, for running model. Do not be alarmed by sizes.')
else
    warning('Emmanuel New Fish 1 data with genes and coords not downloaded on this computer')
    % load('C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/ZebrafishData/FishDataEmmanuel123/Fish1/Fish1cfos.mat')
    % load('C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/ZebrafishData/FishDataEmmanuel123/Fish1/save_Emmanuel_NewFish1_dLDS.mat')
end
%%
plotOption = input('Plot maps (1), plot a gene map (2), plot D column maps (3), or plot histogram of L1/L2 of D columns (4)?');
if plotOption == 1
    percentileThreshold = input('Which quantile? 0.99999 means 20k connections:');
    mapFishConnectivity([],filtered.coords,Phi,F,percentileThreshold)
elseif plotOption == 2
    mapNewFishGenes(filtered.cfos,filtered.coords)
elseif plotOption == 3
    mapDCols(Phi,filtered.coords)
elseif plotOption == 4
    l1l2PerDCol = compareDCols(Phi);
elseif plotOption == 5
    corrCfosToD = corrOneGeneToD(filtered.cfos, Phi);
else
    warning('Select an option above.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% correlate genes to D columns
function corrOneGeneToDOut = corrOneGeneToD(geneData, Phi)
corrOneGeneToDOut = zeros(size(Phi,2),1); % number of columns of Phi
for col = 1:size(Phi,2)
    whos
    r = corrcoef(double(geneData),Phi(:,col));
    corrOneGeneToDOut(col) = (r(1,2))^2;
end
end
%% map genes
function mapNewFishGenes(dataToMap,whichFilteredIdx)
genes = double(dataToMap);
genes(genes==0) = NaN;
figure();histogram(genes)
genes(genes < quantile(genes(:),0.9)) = NaN;
disp('Hardcoded threshold')
for whichGene = 1 %:size(dataToMap,2)
    % f1 = figure();
    figure();
    hold on
    whichCoordsNoNan                            = whichFilteredIdx;
    whichCoordsNoNan(isnan(whichFilteredIdx(:,1)),:) = [];
    shp  = alphaShape(whichCoordsNoNan(:,3),whichCoordsNoNan(:,2),whichCoordsNoNan(:,1));
    % shp1 = criticalAlpha(shp, 'one-region');
    % tri = alphaTriangulation(shp);
    % polyshp = boundaryshape(tri);
    plot(shp,'FaceColor', 'blue', 'FaceAlpha', 0.001, 'EdgeAlpha',0.05)
    axis equal

    
    % title(sprintf('Gene %d',whichGene))
    whichGene = 'cfos';
    disp('Hardcoded - cfos')
    title(whichGene)
    % scatter3(whichFilteredIdx(:,3),whichFilteredIdx(:,2),whichFilteredIdx(:,1),genes(:,whichGene)/10)
    scatter3(whichFilteredIdx(:,3),whichFilteredIdx(:,2),whichFilteredIdx(:,1),genes/10)

    hold off

    myDateTime = datetime('now');
    myDateTime.Format = 'yyyyMMMddHHmmss';
    dtToStr = string(myDateTime);
    filename = sprintf('geneMap_Emmanuel_%s_%s.svg',whichGene,dtToStr);
    print('-painters','-dsvg',filename)

    l1l2gene = sum(abs(genes),'omitnan')./sqrt(sum(genes.^2,'omitnan'));

    fprintf('L1/L2 of gene distribution: %0.9f',l1l2gene)
end

end
%% calculate L1 over L2 for D (observation matrix) columns (clusters of neurons) - Goal: compare distribution of points - concentrated vs. global
function l1l2PerDCol = compareDCols(Phi)
% disp(sum(abs(Phi)))
% disp(sqrt(sum(Phi.^2)))
% disp(sum(Phi.^2,'omitnan')) % all ones - normalized?
% disp(Phi.^2) % very small values
l1l2PerDCol = sum(abs(Phi),'omitnan')./sqrt(sum(Phi.^2,'omitnan'));
[~,idxMaxL1L2Col] = max(l1l2PerDCol);
[~,idxMinL1L2Col] = min(l1l2PerDCol);
% plot histogram
figure();
histogram(l1l2PerDCol)
title('L1/L2 norm of each column of D')
ylabel('Count')
xlabel('L1/L2')

myDateTime = datetime('now');
myDateTime.Format = 'yyyyMMMddHHmmss';
dtToStr = string(myDateTime);
filename = sprintf('histL1L2D_Emmanuel_%s.svg',dtToStr);
print('-painters','-dsvg',filename)
close all

figure();
imagesc(Phi)
colorbar
% set(gca, 'ColorScale', 'log'); 
title(sprintf('D, L1/L2 max col %d min col %d', idxMaxL1L2Col,idxMinL1L2Col))
ylabel('Ambient dimension')
xlabel('Latent dimension')

myDateTime = datetime('now');
myDateTime.Format = 'yyyyMMMddHHmmss';
dtToStr = string(myDateTime);
filename = sprintf('DL1L2_Emmanuel_%s.svg',dtToStr);
print('-painters','-dsvg',filename)
close all

figure();
imagesc(Phi(:,[idxMaxL1L2Col, idxMinL1L2Col]))
colorbar
% set(gca, 'ColorScale', 'log'); 
title(sprintf('D, L1/L2 max col %d min col %d', idxMaxL1L2Col,idxMinL1L2Col))
ylabel('Ambient dimension')
xlabel('Latent dimension')

myDateTime = datetime('now');
myDateTime.Format = 'yyyyMMMddHHmmss';
dtToStr = string(myDateTime);
filename = sprintf('DmaxminL1L2_Emmanuel_%s.svg',dtToStr);
print('-painters','-dsvg',filename)
close all
end
%% map D (observation matrix) columns (clusters of neurons)
function mapDCols(Phi,whichFilteredIdx)
absPhi = abs(Phi);
absPhi(absPhi==0) = NaN;
for col = 1:size(absPhi,2)
    % f1 = figure();
    figure();
    hold on
    whichCoordsNoNan                            = whichFilteredIdx;
    whichCoordsNoNan(isnan(whichFilteredIdx(:,1)),:) = [];
    shp  = alphaShape(whichCoordsNoNan(:,3),whichCoordsNoNan(:,2),whichCoordsNoNan(:,1));
    % shp1 = criticalAlpha(shp, 'one-region');
    % tri = alphaTriangulation(shp);
    % polyshp = boundaryshape(tri);
    plot(shp,'FaceColor', 'blue', 'FaceAlpha', 0.001, 'EdgeAlpha',0.05)
    axis equal

    
    title(sprintf('D column %d (abs. val.)',col))
    % scatter3(whichFilteredIdx(:,3),whichFilteredIdx(:,2),whichFilteredIdx(:,1),genes(:,whichGene)/10)
    scatter3(whichFilteredIdx(:,3),whichFilteredIdx(:,2),whichFilteredIdx(:,1),absPhi(:,col))

    hold off

    myDateTime = datetime('now');
    myDateTime.Format = 'yyyyMMMddHHmmss';
    dtToStr = string(myDateTime);
    filename = sprintf('DcolMap_EmmanuelNewFish1_Col%d_%s.svg',col,dtToStr);
    % filename = erase(filename, ".");
    print('-painters','-dsvg',filename)
end

end

%% map connections (DO adjacency matrices)
function varargout = mapFishConnectivity(~,whichFilteredIdx,Phi,...
    F, whichPercentile)

ifPlotFlows = 0;

whichCoordsNoNan                            = whichFilteredIdx;
whichCoordsNoNan(isnan(whichFilteredIdx(:,1)),:) = [];
shp  = alphaShape(whichCoordsNoNan(:,3),whichCoordsNoNan(:,2),whichCoordsNoNan(:,1));


for whichF = 1:size(F,2)
    figure();
    hold on
    
    plot(shp,'FaceColor', 'blue', 'FaceAlpha', 0.001, 'EdgeAlpha',0.01)
    axis equal

    unitnormF           = F{1,whichF}.*1./max(abs(F{1,whichF}),[],'all');
    connectivityDueToFi = Phi*unitnormF*Phi.';
    % figure();histogram(connectivityDueToFi);
    % pause();
    
    % parpool(4)
    % parfor source = 1:size(connectivityDueToFi,1)
    upperQuantileThresholdStrength = quantile(abs(connectivityDueToFi),whichPercentile,"all");
    thresholdToPlot                = upperQuantileThresholdStrength;

    whichConnAboveThr = abs(connectivityDueToFi) > thresholdToPlot;
    connAboveThr = connectivityDueToFi(whichConnAboveThr);
    connAndIdx1 = repmat(whichFilteredIdx(:,3), 1, size(connectivityDueToFi,1)); % target x loc
    connAndIdx1 = connAndIdx1(whichConnAboveThr);
    connAndIdx2 = repmat(whichFilteredIdx(:,2), 1, size(connectivityDueToFi,1)); % target y loc
    connAndIdx2 = connAndIdx2(whichConnAboveThr);
    connAndIdx3 = repmat(whichFilteredIdx(:,1), 1, size(connectivityDueToFi,1)); % target z loc
    connAndIdx3 = connAndIdx3(whichConnAboveThr);
    connAndIdx4 = repmat(whichFilteredIdx(:,3).', size(connectivityDueToFi,1), 1); % source x loc
    connAndIdx4 = connAndIdx4(whichConnAboveThr);
    connAndIdx5 = repmat(whichFilteredIdx(:,2).', size(connectivityDueToFi,1), 1); % source y loc
    connAndIdx5 = connAndIdx5(whichConnAboveThr);
    connAndIdx6 = repmat(whichFilteredIdx(:,1).', size(connectivityDueToFi,1), 1); % source z loc
    connAndIdx6 = connAndIdx6(whichConnAboveThr);


    for whichConnection = 1:size(connAboveThr,1)
        p1 = [connAndIdx1(whichConnection) connAndIdx2(whichConnection) connAndIdx3(whichConnection)];   % Target point
        p2 = [connAndIdx4(whichConnection) connAndIdx5(whichConnection) connAndIdx6(whichConnection)];       % Source Point

        dp = p1-p2;                         % Difference

        if connAboveThr(whichConnection) > 0
            arrowColor = 'b';
        else
            arrowColor = 'r';
        end

        hold on

        % upperQuantileThresholdStrength = quantile(abs(connectivityDueToFi),0.99,"all");
        
        if sum(dp,'all') ~= 0 %& abs(connectivityDueToFi(source,target)) > thresholdToPlot
            fprintf('\n %d out of %d %.2f', whichConnection, size(connAboveThr,1), abs(connAboveThr(whichConnection)))
            scatter3(p1(1),p1(2),p1(3),50*abs(connAboveThr(whichConnection)),arrowColor)
            scatter3(p2(1),p2(2),p2(3),50*abs(connAboveThr(whichConnection)),arrowColor)
            
            % quiver3(p2(1),p2(2),p2(3),dp(1),dp(2),dp(3),...
            %     'LineWidth',abs(connAboveThr(whichConnection)),...
            %     'ShowArrowHead','on',...
            %     'Color',arrowColor)
        end
        % https://www.mathworks.com/matlabcentral/answers/2083223-plotting-lines-with-quadruplets-r-g-b-alpha
    end
    
    
hold off

myDateTime = datetime('now');
myDateTime.Format = 'yyyyMMMddHHmmss';
dtToStr = string(myDateTime);

filename = sprintf('connecMap_Emmanuel_NewFish1_%.9f_%d_%s.svg',whichPercentile,whichF,dtToStr);
filename = erase(filename, ".");
print('-painters','-dsvg',filename)

end




varargout{1} = f1;
end

%%
function traceVsBlock = lookupTraceVsBlock3D(whichFilteredIdx,bs1,bs2)

% same trace used multiple times across locations on the map, goes from 0 
% to 14089 (originally python syntax, so added 1 to go from 1 to 14090
% instead)
cellIDs     = whichFilteredIdx;
remappedIDs = 1:size(whichFilteredIdx,1);

% store which block(s) each remappedID belongs to
blocksize1 = bs1; % divisible into 512 and 560 - hardcoded - FIXME
blocksize2 = bs2; % divisible into 512 and 560 - hardcoded - FIXME
% blocksize3 = 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%