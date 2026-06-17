%% for the TMLR run of b-dLDS with last 40 samples - main example
%
% just use saveFMT_forCLDS_1bhv.mat and cldsoutput_TMLR_last40trials.mat

% thisSim = load('C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code/saveFMT_1bhvgraded_forCLDS_TMLR40_rerun_GOOD.mat');
% 
% %
% dFF3d = reshape(cell2mat(thisSim.dFF).',  size(thisSim.dFF{1},2),size(thisSim.dFF{1},1),length(thisSim.dFF)); % goes columnwise
% bhv3d = reshape(cell2mat(thisSim.behaviorData).', size(thisSim.behaviorData{1},2),size(thisSim.behaviorData{1},1),length(thisSim.behaviorData));
% gts3d = reshape(cell2mat(thisSim.groundTruthStates).', size(thisSim.groundTruthStates{1},2),size(thisSim.groundTruthStates{1},1),length(thisSim.groundTruthStates));
% 
% % Make ground truth
% 
% Fgt2 = cell(numel(thisSim.Fgt{1})+numel(thisSim.Fgt{2})+1,1);
% for ll = 1:numel(Fgt2);   Fgt2{ll} = zeros(thisSim.sizeD,thisSim.sizeD);              end
% for ll = 1:numel(thisSim.Fgt{1}); Fgt2{ll}(1:thisSim.sizeD/2,1:thisSim.sizeD/2) = thisSim.Fgt{1}{ll}; end
% for ll = 1:numel(thisSim.Fgt{2})
%     Fgt2{ll+numel(thisSim.Fgt{1})}(thisSim.sizeD/2+1:thisSim.sizeD,thisSim.sizeD/2+1:thisSim.sizeD) = thisSim.Fgt{2}{ll}; 
% end
% Fgt2{end}                        = 0.001*randn(thisSim.sizeD,thisSim.sizeD);
% 
% %
% fcMat = zeros(8,8,3000,50); % one time trace for each entry in combined Fc matrix, for each sample
% for i = 1:length(thisSim.groundTruthStates) % samples N_ex
%     for j = 1:length(Fgt2)-1 % true operators
%         thisC  = thisSim.groundTruthStates{i}(j,:); % 1 by 3000 time points
%         thisF  = Fgt2{j};
%         % thisFC = repmat(thisC,[8 8 1]);
%         for k = 1:3000 % time points
%             fcMat(:,:,k,i) = fcMat(:,:,k,i) + thisF*thisC(k);
%         end
%     end
% end
% 
% %
% 
% save(fullfile(folderPathOut,sprintf('forCLDS_%s',fileList(thisFile).name)))

%% for the TMLR runs of b-dLDS with only 40 samples 
folderPathIn  = 'C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code/TMLR';
folderPathOut = 'C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code/TMLROut';
fileExtension = '*.mat'; % Change to your desired file extension (e.g., *.mat, *.csv)
fileList = dir(fullfile(folderPathIn, fileExtension));
%
for thisFile = 1:length(fileList)
    filePath = fullfile(folderPathIn, fileList(thisFile).name);
    thisSim = load(filePath);

    %
    dFF3d = reshape(cell2mat(thisSim.dFF).',  size(thisSim.dFF{1},2),size(thisSim.dFF{1},1),length(thisSim.dFF)); % goes columnwise
    bhv3d = reshape(cell2mat(thisSim.behaviorData).', size(thisSim.behaviorData{1},2),size(thisSim.behaviorData{1},1),length(thisSim.behaviorData));
    gts3d = reshape(cell2mat(thisSim.groundTruthStates).', size(thisSim.groundTruthStates{1},2),size(thisSim.groundTruthStates{1},1),length(thisSim.groundTruthStates));

    % Make ground truth

    Fgt2 = cell(numel(thisSim.Fgt{1})+numel(thisSim.Fgt{2})+1,1);
    for ll = 1:numel(Fgt2);   Fgt2{ll} = zeros(thisSim.sizeD,thisSim.sizeD);              end
    for ll = 1:numel(thisSim.Fgt{1}); Fgt2{ll}(1:thisSim.sizeD/2,1:thisSim.sizeD/2) = thisSim.Fgt{1}{ll}; end
    for ll = 1:numel(thisSim.Fgt{2})
        Fgt2{ll+numel(thisSim.Fgt{1})}(thisSim.sizeD/2+1:thisSim.sizeD,thisSim.sizeD/2+1:thisSim.sizeD) = thisSim.Fgt{2}{ll}; 
    end
    Fgt2{end}                        = 0.001*randn(thisSim.sizeD,thisSim.sizeD);

    %
    fcMat = zeros(8,8,3000,50); % one time trace for each entry in combined Fc matrix, for each sample
    for i = 1:length(thisSim.groundTruthStates) % samples N_ex
        for j = 1:length(Fgt2)-1 % true operators
            thisC  = thisSim.groundTruthStates{i}(j,:); % 1 by 3000 time points
            thisF  = Fgt2{j};
            % thisFC = repmat(thisC,[8 8 1]);
            for k = 1:3000 % time points
                fcMat(:,:,k,i) = fcMat(:,:,k,i) + thisF*thisC(k);
            end
        end
    end

    %

    save(fullfile(folderPathOut,sprintf('forCLDS_%s',fileList(thisFile).name)))
end

%% also for plotting for TMLR: 10 bhv version
folderPathIn  = 'C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code/TMLR10bhv';
folderPathOut = 'C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code/TMLROut';
fileExtension = '*.mat'; % Change to your desired file extension (e.g., *.mat, *.csv)
fileList = dir(fullfile(folderPathIn, fileExtension));
%
for thisFile = 1:length(fileList)
    filePath = fullfile(folderPathIn, fileList(thisFile).name);
    thisSim = load(filePath);

    %
    dFF3d = reshape(cell2mat(thisSim.dFF).',  size(thisSim.dFF{1},2),size(thisSim.dFF{1},1),length(thisSim.dFF)); % goes columnwise
    bhv3d = reshape(cell2mat(thisSim.behaviorData).', size(thisSim.behaviorData{1},2),size(thisSim.behaviorData{1},1),length(thisSim.behaviorData));
    gts3d = reshape(cell2mat(thisSim.groundTruthStates).', size(thisSim.groundTruthStates{1},2),size(thisSim.groundTruthStates{1},1),length(thisSim.groundTruthStates));

    % Make ground truth

    Fgt2 = cell(numel(thisSim.Fgt{1})+numel(thisSim.Fgt{2})+1,1);
    for ll = 1:numel(Fgt2);   Fgt2{ll} = zeros(thisSim.sizeD,thisSim.sizeD);              end
    for ll = 1:numel(thisSim.Fgt{1}); Fgt2{ll}(1:thisSim.sizeD/2,1:thisSim.sizeD/2) = thisSim.Fgt{1}{ll}; end
    for ll = 1:numel(thisSim.Fgt{2})
        Fgt2{ll+numel(thisSim.Fgt{1})}(thisSim.sizeD/2+1:thisSim.sizeD,thisSim.sizeD/2+1:thisSim.sizeD) = thisSim.Fgt{2}{ll}; 
    end
    Fgt2{end}                        = 0.001*randn(thisSim.sizeD,thisSim.sizeD);

    %
    fcMat = zeros(8,8,3000,50); % one time trace for each entry in combined Fc matrix, for each sample
    for i = 1:length(thisSim.groundTruthStates) % samples N_ex
        for j = 1:length(Fgt2)-1 % true operators
            thisC  = thisSim.groundTruthStates{i}(j,:); % 1 by 3000 time points
            thisF  = Fgt2{j};
            % thisFC = repmat(thisC,[8 8 1]);
            for k = 1:3000 % time points
                fcMat(:,:,k,i) = fcMat(:,:,k,i) + thisF*thisC(k);
            end
        end
    end

    %

    save(fullfile(folderPathOut,sprintf('forCLDS_%s',fileList(thisFile).name)))
end


%% in the case that we stick with the AISTATS sims with 50 samples train, convert 10 extra samples for testing for CLDS
folderPathIn  = 'C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code/TMLRRegen';
folderPathOut = 'C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code/TMLROutRegen';
fileExtension = '*.mat'; % Change to your desired file extension (e.g., *.mat, *.csv)
fileList = dir(fullfile(folderPathIn, fileExtension));
%
for thisFile = 1:length(fileList)
    filePath = fullfile(folderPathIn, fileList(thisFile).name);
    thisSim = load(filePath);

    %
    dFF3d = reshape(cell2mat(thisSim.dFFExtra).',  size(thisSim.dFFExtra{1},2),size(thisSim.dFFExtra{1},1),length(thisSim.dFFExtra)); % goes columnwise
    bhv3d = reshape(cell2mat(thisSim.behaviorDataExtra).', size(thisSim.behaviorDataExtra{1},2),size(thisSim.behaviorDataExtra{1},1),length(thisSim.behaviorDataExtra));
    gts3d = reshape(cell2mat(thisSim.groundTruthStatesExtra).', size(thisSim.groundTruthStatesExtra{1},2),size(thisSim.groundTruthStatesExtra{1},1),length(thisSim.groundTruthStatesExtra));
    
    % Make ground truth
    
    Fgt2 = cell(numel(thisSim.Fgt{1})+numel(thisSim.Fgt{2})+1,1);
    for ll = 1:numel(Fgt2);   Fgt2{ll} = zeros(thisSim.sizeD,thisSim.sizeD);              end
    for ll = 1:numel(thisSim.Fgt{1}); Fgt2{ll}(1:thisSim.sizeD/2,1:thisSim.sizeD/2) = thisSim.Fgt{1}{ll}; end
    for ll = 1:numel(thisSim.Fgt{2})
        Fgt2{ll+numel(thisSim.Fgt{1})}(thisSim.sizeD/2+1:thisSim.sizeD,thisSim.sizeD/2+1:thisSim.sizeD) = thisSim.Fgt{2}{ll}; 
    end
    Fgt2{end}                        = 0.001*randn(thisSim.sizeD,thisSim.sizeD);
    
    %
    fcMat = zeros(8,8,3000,50); % one time trace for each entry in combined Fc matrix, for each sample
    for i = 1:length(thisSim.groundTruthStates) % samples N_ex
        for j = 1:length(Fgt2)-1 % true operators
            thisC  = thisSim.groundTruthStates{i}(j,:); % 1 by 3000 time points
            thisF  = Fgt2{j};
            % thisFC = repmat(thisC,[8 8 1]);
            for k = 1:3000 % time points
                fcMat(:,:,k,i) = fcMat(:,:,k,i) + thisF*thisC(k);
            end
        end
    end
    
    %
    
    save(fullfile(folderPathOut,sprintf('forCLDS_%s',fileList(thisFile).name)))
end

