addpath(genpath('.'))
%%
clear; close all; clc
%% initialize b-dLDS table, CLDS table (10 seeds x 6 columns of Phi)
meanMSEbdLDS     = zeros(10,6);
meanMSECLDS      = zeros(10,6);
meanMSEbdLDS10   = zeros(10,6);
%% If doing relative MSE, get zeros in denominator. Need to rewrite as flattened - make sure to check dimensions. 
for howManyC = 1:6
    for whichSeed = 1:3%10
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
        any(isnan(trueFC),'all')
        any(isnan(fcMatbdLDS),'all')
        any(isnan(clds),'all')
        for l = 1:8
            for m = 1:8
                for n = 1:50
                    fprintf('\n%d %d %d %d %d\n',howManyC,whichSeed,l,m,n)
                    thistrueFC  = trueFC(l,m,:,n); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
                    thisbdLDSFC = fcMatbdLDS(l,m,:,n); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
                    thisCLDSAt  = clds(l,m,:,n); % in column major order: row 1 col 1, row 2 col 1, row 3 col 1, etc.
       
                    % calculate R^2
                    % corrbdlds   = corrcoef(thisbdLDSFC,thistrueFC);
                    % corrclds    = corrcoef(thisCLDSFC,thistrueFC);
                    % thisr2bdlds = corrbdlds(1,2)^2;
                    % thisr2clds  = corrclds(1,2)^2;
                    %% calculate relative MSE ||A-\hat{A}||2^2/||A||2^2
                    thisMSEbdLDS = sum((thistrueFC-thisbdLDSFC).^2)/length(thistrueFC);
                    % thisMSEbdLDS2 = immse(thistrueFC,thisbdLDSFC); %checked, identical to above
                    thisMSECLDS  = sum((thistrueFC-thisCLDSAt).^2)/length(thistrueFC); %immse(thistrueFC,thisCLDSAt);
                    bslnMSE      = sum((thistrueFC).^2)/length(thistrueFC);

                    thisMSEbdLDS = thisMSEbdLDS/bslnMSE;
                    thisMSECLDS  = thisMSECLDS/bslnMSE;
                    %% store in b-dLDS, CLDS tables (rows: seeds, columns: number of coefficients, i.e., nonzero columns of Phi)
                    meanMSEbdLDS(whichSeed,howManyC) = meanMSEbdLDS(whichSeed,howManyC) + thisMSEbdLDS/(8*8*50);
                    meanMSECLDS(whichSeed,howManyC)  = meanMSECLDS(whichSeed,howManyC)  + thisMSECLDS/(8*8*50);
                end
            end
        end
    end
end




%% plot figure: one line b-dLDS, one CLDS (plot errorbars - mean, stdev) - x axis: number of coefficients tied to behavior
fig1 = figure();


title('Relative MSE(A) vs. no. of true dynamics tied to behavior')
%mean of each column (across seeds), standard error (sample size = #
%seeds = 10? or seeds*rows*cols*samples = 10*8*8*50 = 32000)
hold on
errorbar(1:6,mean(meanMSEbdLDS),std(meanMSEbdLDS)/sqrt(32000),'DisplayName','b-dLDS, 1 bhv')
errorbar(1:6,mean(meanMSECLDS),std(meanMSECLDS)/sqrt(32000),'DisplayName','CLDS, 1 bhv')
yscale("log")
ylabel('Relative MSE of dynamics matrix A over time')
xlabel('Number of nonzero columns of true Phi (dynamics tied to 1 behavior)')
xlim([0 6.5])
legend
hold off
grid off
box off