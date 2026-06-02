addpath(genpath('.'))
%%
clear; close all; clc
%% initialize b-dLDS table, CLDS table (10 seeds x 6 columns of Phi)
meanMSEbdLDS     = zeros(10,6);
meanMSECLDS      = zeros(10,6);
meanMSEbdLDS10   = zeros(10,6);
%% If doing relative MSE, get zeros in denominator. Need to rewrite as flattened - make sure to check dimensions. 
for howManyC = 1:6
    for whichSeed = 1:10
        fprintf('\n%d %d\n',howManyC,whichSeed)
        load(sprintf('C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code/AISTATS/saveFMT_251124_seed%02d_AISTATS_b1_c%d.mat',whichSeed,howManyC))
        load(sprintf('C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code/AISTATSOut/forCLDS_saveFMT_251124_seed%02d_AISTATS_b1_c%d.mat',whichSeed,howManyC))
        load(sprintf("C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/clds/notebooks/251124_seed%02d_AISTATS_b1_c%d.mat",whichSeed,howManyC))
        %% construct 8x8x3000(x50) for learned dLDS model (lin. comb. f and c)
        fcMatbdLDS = zeros(8,8,3000,50); % one time trace for each entry in combined Fc matrix, for each sample
        for i = 1:length(B_cell) % samples N_ex
            for j = 1:length(F) % learned operators
                thisC  = B_cell{i}(j,:); % 1 by 3000 time points
                thisF  = F{j};
                % thisFC = repmat(thisC,[8 8 1]);
                for k = 1:3000 % time points
                    fcMatbdLDS(:,:,k,i) = fcMatbdLDS(:,:,k,i) + thisF*thisC(k);
                end
            end
        end
        %%
        % cldsPermute  = permute(clds,[2 3 1]);
        % any(isnan(trueFC),'all')
        % any(isnan(fcMatbdLDS),'all')
        % any(isnan(clds),'all')
        
        thistrueFC   = reshape(trueFC,1,[]);
        thisbdLDSFC  = reshape(fcMatbdLDS,1,[]);
        thisCLDSAt   = reshape(clds,1,[]);

        %% calculate relative MSE ||A-\hat{A}||_2^2/||A||_2^2
        thisMSEbdLDS = sum((thistrueFC-thisbdLDSFC).^2)/length(thistrueFC);
        % thisMSEbdLDS2 = immse(thistrueFC,thisbdLDSFC); %checked, identical to above
        thisMSECLDS  = sum((thistrueFC-thisCLDSAt).^2)/length(thistrueFC); %immse(thistrueFC,thisCLDSAt);
        bslnMSE      = sum((thistrueFC).^2)/length(thistrueFC);

        thisMSEbdLDS = thisMSEbdLDS/bslnMSE;
        thisMSECLDS  = thisMSECLDS/bslnMSE;
        %% store in b-dLDS, CLDS tables (rows: seeds, columns: number of coefficients, i.e., nonzero columns of Phi)
        meanMSEbdLDS(whichSeed,howManyC) = thisMSEbdLDS;
        meanMSECLDS(whichSeed,howManyC)  = thisMSECLDS;

        % for l = 1:8
        %     for m = 1:8
        %         for n = 1:50
        %             fprintf('\n%d %d %d %d %d\n',howManyC,whichSeed,l,m,n)
        %             thistrueFC  = trueFC(l,m,:,n); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
        %             thisbdLDSFC = fcMatbdLDS(l,m,:,n); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
        %             thisCLDSAt  = clds(l,m,:,n); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
        % 
        %             % calculate R^2
        %             % corrbdlds   = corrcoef(thisbdLDSFC,thistrueFC);
        %             % corrclds    = corrcoef(thisCLDSFC,thistrueFC);
        %             % thisr2bdlds = corrbdlds(1,2)^2;
        %             % thisr2clds  = corrclds(1,2)^2;
        %             %% calculate relative MSE ||A-\hat{A}||2^2/||A||2^2
        %             thisMSEbdLDS = sum((thistrueFC-thisbdLDSFC).^2)/length(thistrueFC);
        %             % thisMSEbdLDS2 = immse(thistrueFC,thisbdLDSFC); %checked, identical to above
        %             thisMSECLDS  = sum((thistrueFC-thisCLDSAt).^2)/length(thistrueFC); %immse(thistrueFC,thisCLDSAt);
        %             bslnMSE      = sum((thistrueFC).^2)/length(thistrueFC);
        % 
        %             thisMSEbdLDS = thisMSEbdLDS/bslnMSE;
        %             thisMSECLDS  = thisMSECLDS/bslnMSE;
        %             %% store in b-dLDS, CLDS tables (rows: seeds, columns: number of coefficients, i.e., nonzero columns of Phi)
        %             meanMSEbdLDS(whichSeed,howManyC) = meanMSEbdLDS(whichSeed,howManyC) + thisMSEbdLDS/(8*8*50);
        %             meanMSECLDS(whichSeed,howManyC)  = meanMSECLDS(whichSeed,howManyC)  + thisMSECLDS/(8*8*50);
        %         end
        %     end
        % end
    end
end

%% 10 behaviors
for howManyC = 1:6
    for whichSeed = 1:10
        fprintf('\n%d %d\n',howManyC,whichSeed)
        load(sprintf('C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code/AISTATS10bhv/saveFMT_251126_seed%02d_AISTATS_b10_c%d.mat',whichSeed,howManyC))
        load(sprintf('C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/backup_dLDSDiscreteBehavior/dLDS-Discrete-Matlab-Model_FishBehavior/code/AISTATSOut/forCLDS_saveFMT_251126_seed%02d_AISTATS_b10_c%d.mat',whichSeed,howManyC))
        % load(sprintf("C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/clds/notebooks/251126_seed%02d_AISTATS_b10_c%d.mat",whichSeed,howManyC))
        
        
        % construct 8x8x3000(x50) for ground truth dLDS model (lin. comb. f and c)
        Fgt2 = cell(numel(Fgt{1})+numel(Fgt{2})+1,1);
        for ll = 1:numel(Fgt2);   Fgt2{ll} = zeros(sizeD,sizeD);              end
        for ll = 1:numel(Fgt{1}); Fgt2{ll}(1:sizeD/2,1:sizeD/2) = Fgt{1}{ll}; end
        for ll = 1:numel(Fgt{2})
            Fgt2{ll+numel(Fgt{1})}(sizeD/2+1:sizeD,sizeD/2+1:sizeD) = Fgt{2}{ll}; 
        end
        Fgt2{end}                        = 0.001*randn(sizeD,sizeD);
        
        %%
        trueFC = zeros(8,8,3000,50); % one time trace for each entry in combined Fc matrix, for each sample
        for i = 1:length(groundTruthStates) % samples N_ex
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
        fcMatbdLDS = zeros(8,8,3000,50); % one time trace for each entry in combined Fc matrix, for each sample
        for i = 1:length(B_cell) % samples N_ex
            for j = 1:length(F) % learned operators
                thisC  = B_cell{i}(j,:); % 1 by 3000 time points
                thisF  = F{j};
                % thisFC = repmat(thisC,[8 8 1]);
                for k = 1:3000 % time points
                    fcMatbdLDS(:,:,k,i) = fcMatbdLDS(:,:,k,i) + thisF*thisC(k);
                end
            end
        end
        %
        % cldsPermute  = permute(clds,[2 3 1]);
        any(isnan(trueFC),'all')
        any(isnan(fcMatbdLDS),'all')
        any(isnan(clds),'all')


        thistrueFC   = reshape(trueFC,1,[]);
        thisbdLDSFC  = reshape(fcMatbdLDS,1,[]);
        % thisCLDSAt   = reshape(clds,1,[]);

        %% calculate relative MSE ||A-\hat{A}||_2^2/||A||_2^2
        thisMSEbdLDS10 = sum((thistrueFC-thisbdLDSFC).^2)/length(thistrueFC);
        % thisMSEbdLDS2 = immse(thistrueFC,thisbdLDSFC); %checked, identical to above
        % thisMSECLDS10  = sum((thistrueFC-thisCLDSAt).^2)/length(thistrueFC); %immse(thistrueFC,thisCLDSAt);
        bslnMSE10      = sum((thistrueFC).^2)/length(thistrueFC);

        thisMSEbdLDS10 = thisMSEbdLDS10/bslnMSE10;
        % thisMSECLDS10  = thisMSECLDS10/bslnMSE10;
        %% store in b-dLDS, CLDS tables (rows: seeds, columns: number of coefficients, i.e., nonzero columns of Phi)
        meanMSEbdLDS10(whichSeed,howManyC) = thisMSEbdLDS10;
        % meanMSECLDS10(whichSeed,howManyC)  = thisMSECLDS10;

        % for l = 1:8
        %     for m = 1:8
        %         for n = 1:50
        %             fprintf('\n%d %d %d %d %d\n',howManyC,whichSeed,l,m,n)
        %             thistrueFC  = trueFC(l,m,:,n); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
        %             thisbdLDSFC = fcMatbdLDS(l,m,:,n); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
        % 
        %             %% calculate relative MSE
        % 
        %             thisMSEbdLDS   = sum((thistrueFC-thisbdLDSFC).^2)/length(thistrueFC);
        %             bslnMSE        = sum((thistrueFC).^2)/length(thistrueFC);
        % 
        %             thisMSEbdLDS10 = thisMSEbdLDS/bslnMSE;
        %             % store in b-dLDS, CLDS tables (rows: seeds, columns: number of coefficients, i.e., nonzero columns of Phi)
        %             meanMSEbdLDS10(whichSeed,howManyC) = meanMSEbdLDS10(whichSeed,howManyC) + thisMSEbdLDS10/(8*8*50);
        %         end
        %     end
        % end
    end
end


%% plot figure: one line b-dLDS, one CLDS (plot errorbars - mean, stdev) - x axis: number of coefficients tied to behavior
fig1 = figure();


%title('Relative MSE(A) vs. no. of true dynamics tied to behavior')
%mean of each column (across seeds), standard error (sample size = #
%seeds = 10? or seeds*rows*cols*samples = 10*8*8*50 = 32000
%Now rMSE calculated across whole table, per model and seed --> averaging across seeds)
hold on
errorbar(1:6,mean(meanMSEbdLDS),std(meanMSEbdLDS)/sqrt(10),'DisplayName','b-dLDS, 1 bhv','LineStyle', 'none', 'Marker', 'o')
errorbar(1:6,mean(meanMSEbdLDS10),std(meanMSEbdLDS10)/sqrt(10),'DisplayName','b-dLDS, 10 bhv','LineStyle', 'none', 'Marker', 'o')
errorbar(1:6,mean(meanMSECLDS),std(meanMSECLDS)/sqrt(10),'DisplayName','CLDS, 1 bhv','LineStyle', 'none', 'Marker', 'o')
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