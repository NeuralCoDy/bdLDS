addpath(genpath('.'))
addpath(genpath('./External_Packages/distinguishable_colors/'))
%%
clear; close all; clc
%% initialize b-dLDS table, CLDS table (10 seeds x 6 columns of Phi)
meanMSEbdLDS      = zeros(10,6);
meanMSEAblFrob    = zeros(10,6);
meanMSEAblSpr     = zeros(10,6);
meanMSEAblBoth    = zeros(10,6);


% meanMSECLDS      = zeros(10,6);
meanMSEbdLDS10   = zeros(10,6);
meanMSEAblFrob10 = zeros(10,6);
meanMSEAblSpr10  = zeros(10,6);
meanMSEAblBoth10 = zeros(10,6);
%% If doing relative MSE, get zeros in denominator. Need to rewrite as flattened - make sure to check dimensions.
% 40/10 train/test split
for howManyC = 1:6
    for whichSeed = 1:10
        fprintf('\n%d %d\n',howManyC,whichSeed)
        % disp('FIXME: check variable names in each file')
        load(sprintf('/cis/home/eyezere1/my_documents/dLDS_fromlaptop_250715/TMLR/saveFMT_260303_seed%02d_TMLR_b1_c%d.mat',whichSeed,howManyC));
        bdLDSAblFrob = load(sprintf('/cis/home/eyezere1/my_documents/dLDS_fromlaptop_250715/Ablation/saveFMT_260407_seed%02d_AblFrobPsi_b1_c%d.mat',whichSeed,howManyC)); % ablate the Frobenius norm on Psi
        bdLDSAblSpr  = load(sprintf('/cis/home/eyezere1/my_documents/dLDS_fromlaptop_250715/Ablation/saveFMT_260407_seed%02d_AblSprPsi_b1_c%d.mat',whichSeed,howManyC)); % ablate the "sparse suspected" setting on Psi, which rescales the step size
        bdLDSAblBoth = load(sprintf('/cis/home/eyezere1/my_documents/dLDS_fromlaptop_250715/Ablation/saveFMT_260416_seed%02d_AblBothFrSp_b1_c%d.mat',whichSeed,howManyC)); % ablate both
        % load(sprintf('C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code/TMLROut/forCLDS_saveFMT_260303_seed%02d_TMLR_b1_c%d.mat',whichSeed,howManyC))
        % load(sprintf("C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/clds/notebooks/260303_seed%02d_TMLR_b1_c%d.mat",whichSeed,howManyC)) 
        
        % construct 8x8x3000(x50) for ground truth dLDS model (lin. comb. f and c)
        Fgt2 = cell(numel(Fgt{1})+numel(Fgt{2})+1,1);
        for ll = 1:numel(Fgt2);   Fgt2{ll} = zeros(sizeD,sizeD);              end
        for ll = 1:numel(Fgt{1}); Fgt2{ll}(1:sizeD/2,1:sizeD/2) = Fgt{1}{ll}; end
        for ll = 1:numel(Fgt{2})
            Fgt2{ll+numel(Fgt{1})}(sizeD/2+1:sizeD,sizeD/2+1:sizeD) = Fgt{2}{ll}; 
        end
        Fgt2{end}                        = 0.001*randn(sizeD,sizeD);
        
        %
        trueFC = zeros(8,8,3000,40); % one time trace for each entry in combined Fc matrix, for each sample
        for i = 1:40%length(groundTruthStates) % samples N_ex
            for j = 1:length(Fgt2)-1 % true operators
                thisC  = groundTruthStates{i}(j,:); % 1 by 3000 time points
                thisF  = Fgt2{j};
                % thisFC = repmat(thisC,[8 8 1]);
                for k = 1:3000 % time points
                    trueFC(:,:,k,i) = trueFC(:,:,k,i) + thisF*thisC(k);
                end
            end
        end
        
        
        % construct 8x8x3000(x50) for learned dLDS model (lin. comb. f and c)
        fcMatbdLDS = zeros(8,8,3000,40); % one time trace for each entry in combined Fc matrix, for each sample
        for i = 1:40%length(B_cell) % samples N_ex
            for j = 1:length(F) % learned operators
                thisC  = B_cell{i}(j,:); % 1 by 3000 time points
                thisF  = F{j};
                % thisFC = repmat(thisC,[8 8 1]);
                for k = 1:3000 % time points
                    fcMatbdLDS(:,:,k,i) = fcMatbdLDS(:,:,k,i) + thisF*thisC(k);
                end
            end
        end

        


        fcMatAblFrob = zeros(8,8,3000,40); % one time trace for each entry in combined Fc matrix, for each sample
        for i = 1:40%length(B_cell) % samples N_ex
            for j = 1:length(F) % learned operators
                thisC  = bdLDSAblFrob.B_cell{i}(j,:); % 1 by 3000 time points
                thisF  = bdLDSAblFrob.F{j};
                % thisFC = repmat(thisC,[8 8 1]);
                for k = 1:3000 % time points
                    fcMatAblFrob(:,:,k,i) = fcMatAblFrob(:,:,k,i) + thisF*thisC(k);
                end
            end
        end

        


        fcMatAblSpr = zeros(8,8,3000,40); % one time trace for each entry in combined Fc matrix, for each sample
        for i = 1:40%length(B_cell) % samples N_ex
            for j = 1:length(F) % learned operators
                thisC  = bdLDSAblSpr.B_cell{i}(j,:); % 1 by 3000 time points
                thisF  = bdLDSAblSpr.F{j};
                % thisFC = repmat(thisC,[8 8 1]);
                for k = 1:3000 % time points
                    fcMatAblSpr(:,:,k,i) = fcMatAblSpr(:,:,k,i) + thisF*thisC(k);
                end
            end
        end

        fcMatAblBoth = zeros(8,8,3000,40); % one time trace for each entry in combined Fc matrix, for each sample
        for i = 1:40%length(B_cell) % samples N_ex
            for j = 1:length(F) % learned operators
                thisC  = bdLDSAblBoth.B_cell{i}(j,:); % 1 by 3000 time points
                thisF  = bdLDSAblBoth.F{j};
                % thisFC = repmat(thisC,[8 8 1]);
                for k = 1:3000 % time points
                    fcMatAblBoth(:,:,k,i) = fcMatAblBoth(:,:,k,i) + thisF*thisC(k);
                end
            end
        end

        
        %
        % cldsPermute  = permute(clds,[2 3 1]);
        any(isnan(trueFC),'all')
        any(isnan(fcMatAblFrob),'all')
        any(isnan(fcMatAblSpr),'all')
        any(isnan(fcMatAblBoth),'all')
        % any(isnan(clds),'all')
        
        thistrueFC     = reshape(trueFC(:,:,:,1:40),1,[]);
        thisbdLDSFC    = reshape(fcMatbdLDS,1,[]);
        thisAblFrobFC  = reshape(fcMatAblFrob,1,[]);
        thisAblSprFC   = reshape(fcMatAblSpr,1,[]);
        thisAblBothFC  = reshape(fcMatAblBoth,1,[]);
        % thisCLDSAt   = reshape(clds(:,:,:,1:40),1,[]);

        % calculate relative MSE ||A-\hat{A}||_2^2/||A||_2^2
        thisMSEbdLDS   = sum((thistrueFC-thisbdLDSFC).^2);
        thisMSEAblFrob = sum((thistrueFC-thisAblFrobFC).^2);
        thisMSEAblSpr  = sum((thistrueFC-thisAblSprFC).^2);
        thisMSEAblBoth = sum((thistrueFC-thisAblBothFC).^2);

        % thisMSEbdLDS2 = immse(thistrueFC,thisbdLDSFC); %checked, identical to above
        % thisMSECLDS  = sum((thistrueFC-thisCLDSAt).^2); %immse(thistrueFC,thisCLDSAt);
        
        bslnMSE      = sum((thistrueFC).^2);

        thisMSEbdLDS   = thisMSEbdLDS/bslnMSE;
        thisMSEAblFrob = thisMSEAblFrob/bslnMSE;
        thisMSEAblSpr  = thisMSEAblSpr/bslnMSE;
        thisMSEAblBoth = thisMSEAblBoth/bslnMSE;
        % thisMSECLDS  = thisMSECLDS/bslnMSE;
        % store in b-dLDS, CLDS tables (rows: seeds, columns: number of coefficients, i.e., nonzero columns of Phi)
        meanMSEbdLDS(whichSeed,howManyC)   = thisMSEbdLDS;
        meanMSEAblFrob(whichSeed,howManyC) = thisMSEAblFrob;
        meanMSEAblSpr(whichSeed,howManyC)  = thisMSEAblSpr;
        meanMSEAblBoth(whichSeed,howManyC) = thisMSEAblBoth;
        % meanMSECLDS(whichSeed,howManyC)  = thisMSECLDS;

        
    end
end

%% 10 behaviors
for howManyC = 1:6
    for whichSeed = 1:10
        fprintf('\n%d %d\n',howManyC,whichSeed)
        load(sprintf('/cis/home/eyezere1/my_documents/dLDS_fromlaptop_250715/TMLR10bhv/saveFMT_260303_seed%02d_TMLR_b10_c%d.mat',whichSeed,howManyC))
        bdLDSAblFrob = load(sprintf('/cis/home/eyezere1/my_documents/dLDS_fromlaptop_250715/Ablation/saveFMT_260417_seed%02d_AblFrobPsi_b10_c%d.mat',whichSeed,howManyC)); % ablate the Frobenius norm on Psi
        bdLDSAblSpr  = load(sprintf('/cis/home/eyezere1/my_documents/dLDS_fromlaptop_250715/Ablation/saveFMT_260417_seed%02d_AblSprPsi_b10_c%d.mat',whichSeed,howManyC)); % ablate the "sparse suspected" setting on Psi, which rescales the step size
        bdLDSAblBoth = load(sprintf('/cis/home/eyezere1/my_documents/dLDS_fromlaptop_250715/Ablation/saveFMT_260417_seed%02d_AblBothFrSp_b10_c%d.mat',whichSeed,howManyC)); % ablate both
        % load(sprintf('C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code/TMLROut/forCLDS_saveFMT_260303_seed%02d_TMLR_b10_c%d.mat',whichSeed,howManyC))
        % load(sprintf("C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/clds/notebooks/251126_seed%02d_AISTATS_b10_c%d.mat",whichSeed,howManyC))



        % construct 8x8x3000(x50) for ground truth dLDS model (lin. comb. f and c)
        Fgt2 = cell(numel(Fgt{1})+numel(Fgt{2})+1,1);
        for ll = 1:numel(Fgt2);   Fgt2{ll} = zeros(sizeD,sizeD);              end
        for ll = 1:numel(Fgt{1}); Fgt2{ll}(1:sizeD/2,1:sizeD/2) = Fgt{1}{ll}; end
        for ll = 1:numel(Fgt{2})
            Fgt2{ll+numel(Fgt{1})}(sizeD/2+1:sizeD,sizeD/2+1:sizeD) = Fgt{2}{ll}; 
        end
        Fgt2{end}                        = 0.001*randn(sizeD,sizeD);

        %
        trueFC = zeros(8,8,3000,40); % one time trace for each entry in combined Fc matrix, for each sample
        for i = 1:40%length(groundTruthStates) % samples N_ex
            for j = 1:length(Fgt2)-1 % true operators
                thisC  = groundTruthStates{i}(j,:); % 1 by 3000 time points
                thisF  = Fgt2{j};
                % thisFC = repmat(thisC,[8 8 1]);
                for k = 1:3000 % time points
                    trueFC(:,:,k,i) = trueFC(:,:,k,i) + thisF*thisC(k);
                end
            end
        end

        % construct 8x8x3000(x50) for learned dLDS model (lin. comb. f and c)
        fcMatbdLDS = zeros(8,8,3000,40); % one time trace for each entry in combined Fc matrix, for each sample
        for i = 1:40%length(B_cell) % samples N_ex
            for j = 1:length(F) % learned operators
                thisC  = B_cell{i}(j,:); % 1 by 3000 time points
                thisF  = F{j};
                % thisFC = repmat(thisC,[8 8 1]);
                for k = 1:3000 % time points
                    fcMatbdLDS(:,:,k,i) = fcMatbdLDS(:,:,k,i) + thisF*thisC(k);
                end
            end
        end

        fcMatAblFrob = zeros(8,8,3000,40); % one time trace for each entry in combined Fc matrix, for each sample
        for i = 1:40%length(B_cell) % samples N_ex
            for j = 1:length(F) % learned operators
                thisC  = bdLDSAblFrob.B_cell{i}(j,:); % 1 by 3000 time points
                thisF  = bdLDSAblFrob.F{j};
                % thisFC = repmat(thisC,[8 8 1]);
                for k = 1:3000 % time points
                    fcMatAblFrob(:,:,k,i) = fcMatAblFrob(:,:,k,i) + thisF*thisC(k);
                end
            end
        end

        


        fcMatAblSpr = zeros(8,8,3000,40); % one time trace for each entry in combined Fc matrix, for each sample
        for i = 1:40%length(B_cell) % samples N_ex
            for j = 1:length(F) % learned operators
                thisC  = bdLDSAblSpr.B_cell{i}(j,:); % 1 by 3000 time points
                thisF  = bdLDSAblSpr.F{j};
                % thisFC = repmat(thisC,[8 8 1]);
                for k = 1:3000 % time points
                    fcMatAblSpr(:,:,k,i) = fcMatAblSpr(:,:,k,i) + thisF*thisC(k);
                end
            end
        end

        fcMatAblBoth = zeros(8,8,3000,40); % one time trace for each entry in combined Fc matrix, for each sample
        for i = 1:40%length(B_cell) % samples N_ex
            for j = 1:length(F) % learned operators
                thisC  = bdLDSAblBoth.B_cell{i}(j,:); % 1 by 3000 time points
                thisF  = bdLDSAblBoth.F{j};
                % thisFC = repmat(thisC,[8 8 1]);
                for k = 1:3000 % time points
                    fcMatAblBoth(:,:,k,i) = fcMatAblBoth(:,:,k,i) + thisF*thisC(k);
                end
            end
        end



        any(isnan(fcMatAblFrob),'all')
        any(isnan(fcMatAblSpr),'all')
        any(isnan(fcMatAblBoth),'all')
        % any(isnan(clds),'all')
        
        thistrueFC     = reshape(trueFC(:,:,:,1:40),1,[]);
        thisbdLDSFC    = reshape(fcMatbdLDS,1,[]);
        thisAblFrobFC  = reshape(fcMatAblFrob,1,[]);
        thisAblSprFC   = reshape(fcMatAblSpr,1,[]);
        thisAblBothFC  = reshape(fcMatAblBoth,1,[]);
        % thisCLDSAt   = reshape(clds(:,:,:,1:40),1,[]);

        % calculate relative MSE ||A-\hat{A}||_F^2/||A||_2^2
        thisMSEbdLDS   = sum((thistrueFC-thisbdLDSFC).^2);
        thisMSEAblFrob = sum((thistrueFC-thisAblFrobFC).^2);
        thisMSEAblSpr  = sum((thistrueFC-thisAblSprFC).^2);
        thisMSEAblBoth = sum((thistrueFC-thisAblBothFC).^2);

        % thisMSEbdLDS2 = immse(thistrueFC,thisbdLDSFC); %checked, identical to above
        % thisMSECLDS  = sum((thistrueFC-thisCLDSAt).^2); %immse(thistrueFC,thisCLDSAt);
        
        bslnMSE      = sum((thistrueFC).^2);

        thisMSEbdLDS   = thisMSEbdLDS/bslnMSE;
        thisMSEAblFrob = thisMSEAblFrob/bslnMSE;
        thisMSEAblSpr  = thisMSEAblSpr/bslnMSE;
        thisMSEAblBoth = thisMSEAblBoth/bslnMSE;
        % thisMSECLDS  = thisMSECLDS/bslnMSE;
        % store in b-dLDS, CLDS tables (rows: seeds, columns: number of coefficients, i.e., nonzero columns of Phi)
        meanMSEbdLDS10(whichSeed,howManyC)   = thisMSEbdLDS;
        meanMSEAblFrob10(whichSeed,howManyC) = thisMSEAblFrob;
        meanMSEAblSpr10(whichSeed,howManyC)  = thisMSEAblSpr;
        meanMSEAblBoth10(whichSeed,howManyC) = thisMSEAblBoth;
        % meanMSECLDS(whichSeed,howManyC)  = thisMSECLDS;


    end
end


%% plot figure: one line b-dLDS, one CLDS (plot errorbars - mean, stdev) - x axis: number of coefficients tied to behavior
fig1 = figure();


%title('Relative MSE(A) vs. no. of true dynamics tied to behavior')
%mean of each column (across seeds), standard error (sample size = #
%seeds = 10? or seeds*rows*cols*samples = 10*8*8*50 = 32000
%Now rMSE calculated across whole table, per model and seed --> averaging across seeds)
hold on

cmapMaxDistinct = distinguishable_colors(8);

errorbar(1:6,mean(meanMSEbdLDS),std(meanMSEbdLDS)/sqrt(10),'DisplayName','b-dLDS, 1 bhv','LineStyle', 'none', 'Marker', 'o','Color',cmapMaxDistinct(1,:))
errorbar(1:6,mean(meanMSEAblFrob),std(meanMSEAblFrob)/sqrt(10),'DisplayName','Ablate Frobenius norm on Psi, 1 bhv','LineStyle', 'none', 'Marker', 'o','Color',cmapMaxDistinct(2,:))
errorbar(1:6,mean(meanMSEAblSpr),std(meanMSEAblSpr)/sqrt(10),'DisplayName','Ablate Psi step size rescale, 1 bhv','LineStyle', 'none', 'Marker', 'o','Color',cmapMaxDistinct(3,:))
errorbar(1:6,mean(meanMSEAblBoth),std(meanMSEAblBoth)/sqrt(10),'DisplayName','Ablate both, 1 bhv','LineStyle', 'none', 'Marker', 'o','Color',cmapMaxDistinct(4,:))

errorbar(1:6,mean(meanMSEbdLDS10),std(meanMSEbdLDS10)/sqrt(10),'DisplayName','b-dLDS, 10 bhv','LineStyle', 'none', 'Marker', 'o','Color',cmapMaxDistinct(5,:))
errorbar(1:6,mean(meanMSEAblFrob10),std(meanMSEAblFrob10)/sqrt(10),'DisplayName','Ablate Frobenius norm on Psi, 10 bhv','LineStyle', 'none', 'Marker', 'o','Color',cmapMaxDistinct(6,:))
errorbar(1:6,mean(meanMSEAblSpr10),std(meanMSEAblSpr10)/sqrt(10),'DisplayName','Ablate Psi step size rescale, 10 bhv','LineStyle', 'none', 'Marker', 'o','Color',cmapMaxDistinct(7,:))
errorbar(1:6,mean(meanMSEAblBoth10),std(meanMSEAblBoth10)/sqrt(10),'DisplayName','Ablate both, 10 bhv','LineStyle', 'none', 'Marker', 'o','Color',cmapMaxDistinct(8,:))



% errorbar(1:6,mean(meanMSEbdLDS10),std(meanMSEbdLDS10)/sqrt(10),'DisplayName','b-dLDS, 10 bhv','LineStyle', 'none', 'Marker', 'o')
% errorbar(1:6,mean(meanMSECLDS),std(meanMSECLDS)/sqrt(10),'DisplayName','CLDS, 1 bhv','LineStyle', 'none', 'Marker', 'o')
% errorbar(1:6,mean(meanMSECLDS10),std(meanMSECLDS10)/sqrt(10),'DisplayName','CLDS, 10 bhv','LineStyle', 'none', 'Marker', 'o')

% yscale("log")
ylabel('Relative MSE of dynamics matrix A')
xlabel('No. of dynamics tied to behavior)')
xlim([0 6.5])
% xticklabels({sprintf('%0.1f%%',100*1/6),sprintf('%0.1f%%',100*2/6),sprintf('%0.1f%%',100*3/6),sprintf('%0.1f%%',100*4/6),sprintf('%0.1f%%',100*5/6),sprintf('%0.1f%%',100*6/6)})
legend
hold off
grid off
box off

myDateTime = datetime('now');
myDateTime.Format = 'yyyyMMMddHHmmss';
dtToStr = string(myDateTime);

filenamerMSEfig = sprintf('AblationComp_%s.fig',dtToStr);
filenamerMSE    = sprintf('AblationComp_%s.mat',dtToStr);


savefig(fig1, filenamerMSEfig)
save(filenamerMSE, 'meanMSEbdLDS', 'meanMSEAblFrob','meanMSEAblSpr','meanMSEAblBoth')