%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

% parpool(16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code
addpath(genpath('.'))
%% Overall Phi dimensions for this 10-seed test
nBhv = input('How many behavior traces? >=1:');
howManyColsPhi = input('How many columns of Phi (dyn coefs) related to behavior? <=6:');
todaysDate = input('Date (e.g., 251126):');

%%
rng('shuffle')
currentState01 = rng;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate data
% use 'sizeD' latent states, 'sizeD' "recorded channels", 'sizeT' timepoints
load(sprintf('./AISTATS/saveFMT_251124_seed01_AISTATS_b%d_c%d.mat',nBhv,howManyColsPhi))

[latentStatesXExtra, FgtNew, groundTruthStatesExtra] = generateContMultiSubNetwork_Range_Regenerate(...
               sizeT,[sizeD/2,sizeD/2],nF,50,contOpt,false,false,true,Fgt);

%%
% nBhv = 1;
A = randn(nBhv,sum(nF,'all')); 
if howManyColsPhi < 6
    A(:,howManyColsPhi+1:end) = 0; 
end
simulatedPsi = A * 1./vecnorm(A);
simulatedPsi(isnan(simulatedPsi))=0;

justOneGTStateExtra = groundTruthStatesExtra;
behaviorDataExtra = cellfun(@(x) simulatedPsi*x,justOneGTState,'UniformOutput',false);
%%
dFFExtra = statesToObservations(simulatedDmatrix, latentStatesX, noiseLevel);
save(sprintf('./TMLR/saveFMT_%d_seed01_TMLR_b%d_c%d_regen.mat',todaysDate,nBhv,howManyColsPhi))
%%
rng('shuffle')
currentState02 = rng;
%%
load(sprintf('./AISTATS/saveFMT_251124_seed02_AISTATS_b%d_c%d.mat',nBhv,howManyColsPhi))

[latentStatesXExtra, FgtNew, groundTruthStatesExtra] = generateContMultiSubNetwork_Range_Regenerate(...
               sizeT,[sizeD/2,sizeD/2],nF,50,contOpt,false,false,true,Fgt);

%%
% nBhv = 1;
A = randn(nBhv,sum(nF,'all')); 
if howManyColsPhi < 6
    A(:,howManyColsPhi+1:end) = 0; 
end
simulatedPsi = A * 1./vecnorm(A);
simulatedPsi(isnan(simulatedPsi))=0;

justOneGTStateExtra = groundTruthStatesExtra;
behaviorDataExtra = cellfun(@(x) simulatedPsi*x,justOneGTState,'UniformOutput',false);
%%
dFFExtra = statesToObservations(simulatedDmatrix, latentStatesX, noiseLevel);
save(sprintf('./TMLR/saveFMT_%d_seed02_TMLR_b%d_c%d_regen.mat',todaysDate,nBhv,howManyColsPhi))
%%
rng('shuffle')
currentState03 = rng;
%%
load(sprintf('./AISTATS/saveFMT_251124_seed03_AISTATS_b%d_c%d.mat',nBhv,howManyColsPhi))

[latentStatesXExtra, FgtNew, groundTruthStatesExtra] = generateContMultiSubNetwork_Range_Regenerate(...
               sizeT,[sizeD/2,sizeD/2],nF,50,contOpt,false,false,true,Fgt);

%%
% nBhv = 1;
A = randn(nBhv,sum(nF,'all')); 
if howManyColsPhi < 6
    A(:,howManyColsPhi+1:end) = 0; 
end
simulatedPsi = A * 1./vecnorm(A);
simulatedPsi(isnan(simulatedPsi))=0;

justOneGTStateExtra = groundTruthStatesExtra;
behaviorDataExtra = cellfun(@(x) simulatedPsi*x,justOneGTState,'UniformOutput',false);
%%
dFFExtra = statesToObservations(simulatedDmatrix, latentStatesX, noiseLevel);
save(sprintf('./TMLR/saveFMT_%d_seed03_TMLR_b%d_c%d_regen.mat',todaysDate,nBhv,howManyColsPhi))
%%
rng('shuffle')
currentState04 = rng;
%%
load(sprintf('./AISTATS/saveFMT_251124_seed04_AISTATS_b%d_c%d.mat',nBhv,howManyColsPhi))

[latentStatesXExtra, FgtNew, groundTruthStatesExtra] = generateContMultiSubNetwork_Range_Regenerate(...
               sizeT,[sizeD/2,sizeD/2],nF,50,contOpt,false,false,true,Fgt);

%%
% nBhv = 1;
A = randn(nBhv,sum(nF,'all')); 
if howManyColsPhi < 6
    A(:,howManyColsPhi+1:end) = 0; 
end
simulatedPsi = A * 1./vecnorm(A);
simulatedPsi(isnan(simulatedPsi))=0;

justOneGTStateExtra = groundTruthStatesExtra;
behaviorDataExtra = cellfun(@(x) simulatedPsi*x,justOneGTState,'UniformOutput',false);
%%
dFFExtra = statesToObservations(simulatedDmatrix, latentStatesX, noiseLevel);
save(sprintf('./TMLR/saveFMT_%d_seed04_TMLR_b%d_c%d_regen.mat',todaysDate,nBhv,howManyColsPhi))
%%
rng('shuffle')
currentState05 = rng;
%%
load(sprintf('./AISTATS/saveFMT_251124_seed05_AISTATS_b%d_c%d.mat',nBhv,howManyColsPhi))

[latentStatesXExtra, FgtNew, groundTruthStatesExtra] = generateContMultiSubNetwork_Range_Regenerate(...
               sizeT,[sizeD/2,sizeD/2],nF,50,contOpt,false,false,true,Fgt);

%%
% nBhv = 1;
A = randn(nBhv,sum(nF,'all')); 
if howManyColsPhi < 6
    A(:,howManyColsPhi+1:end) = 0; 
end
simulatedPsi = A * 1./vecnorm(A);
simulatedPsi(isnan(simulatedPsi))=0;

justOneGTStateExtra = groundTruthStatesExtra;
behaviorDataExtra = cellfun(@(x) simulatedPsi*x,justOneGTState,'UniformOutput',false);
%%
dFFExtra = statesToObservations(simulatedDmatrix, latentStatesX, noiseLevel);
save(sprintf('./TMLR/saveFMT_%d_seed05_TMLR_b%d_c%d_regen.mat',todaysDate,nBhv,howManyColsPhi))
%%
rng('shuffle')
currentState06 = rng;
%%
load(sprintf('./AISTATS/saveFMT_251124_seed06_AISTATS_b%d_c%d.mat',nBhv,howManyColsPhi))

[latentStatesXExtra, FgtNew, groundTruthStatesExtra] = generateContMultiSubNetwork_Range_Regenerate(...
               sizeT,[sizeD/2,sizeD/2],nF,50,contOpt,false,false,true,Fgt);

%%
% nBhv = 1;
A = randn(nBhv,sum(nF,'all')); 
if howManyColsPhi < 6
    A(:,howManyColsPhi+1:end) = 0; 
end
simulatedPsi = A * 1./vecnorm(A);
simulatedPsi(isnan(simulatedPsi))=0;

justOneGTStateExtra = groundTruthStatesExtra;
behaviorDataExtra = cellfun(@(x) simulatedPsi*x,justOneGTState,'UniformOutput',false);
%%
dFFExtra = statesToObservations(simulatedDmatrix, latentStatesX, noiseLevel);
save(sprintf('./TMLR/saveFMT_%d_seed06_TMLR_b%d_c%d_regen.mat',todaysDate,nBhv,howManyColsPhi))
%%
rng('shuffle')
currentState07 = rng;
%%
load(sprintf('./AISTATS/saveFMT_251124_seed07_AISTATS_b%d_c%d.mat',nBhv,howManyColsPhi))

[latentStatesXExtra, FgtNew, groundTruthStatesExtra] = generateContMultiSubNetwork_Range_Regenerate(...
               sizeT,[sizeD/2,sizeD/2],nF,50,contOpt,false,false,true,Fgt);

%%
% nBhv = 1;
A = randn(nBhv,sum(nF,'all')); 
if howManyColsPhi < 6
    A(:,howManyColsPhi+1:end) = 0; 
end
simulatedPsi = A * 1./vecnorm(A);
simulatedPsi(isnan(simulatedPsi))=0;

justOneGTStateExtra = groundTruthStatesExtra;
behaviorDataExtra = cellfun(@(x) simulatedPsi*x,justOneGTState,'UniformOutput',false);
%%
dFFExtra = statesToObservations(simulatedDmatrix, latentStatesX, noiseLevel);
save(sprintf('./TMLR/saveFMT_%d_seed07_TMLR_b%d_c%d_regen.mat',todaysDate,nBhv,howManyColsPhi))
%%
rng('shuffle')
currentState08 = rng;
%%
load(sprintf('./AISTATS/saveFMT_251124_seed08_AISTATS_b%d_c%d.mat',nBhv,howManyColsPhi))

[latentStatesXExtra, FgtNew, groundTruthStatesExtra] = generateContMultiSubNetwork_Range_Regenerate(...
               sizeT,[sizeD/2,sizeD/2],nF,50,contOpt,false,false,true,Fgt);

%%
% nBhv = 1;
A = randn(nBhv,sum(nF,'all')); 
if howManyColsPhi < 6
    A(:,howManyColsPhi+1:end) = 0; 
end
simulatedPsi = A * 1./vecnorm(A);
simulatedPsi(isnan(simulatedPsi))=0;

justOneGTStateExtra = groundTruthStatesExtra;
behaviorDataExtra = cellfun(@(x) simulatedPsi*x,justOneGTState,'UniformOutput',false);
%%
dFFExtra = statesToObservations(simulatedDmatrix, latentStatesX, noiseLevel);
save(sprintf('./TMLR/saveFMT_%d_seed08_TMLR_b%d_c%d_regen.mat',todaysDate,nBhv,howManyColsPhi))
%%
rng('shuffle')
currentState09 = rng;
%%
load(sprintf('./AISTATS/saveFMT_251124_seed09_AISTATS_b%d_c%d.mat',nBhv,howManyColsPhi))

[latentStatesXExtra, FgtNew, groundTruthStatesExtra] = generateContMultiSubNetwork_Range_Regenerate(...
               sizeT,[sizeD/2,sizeD/2],nF,50,contOpt,false,false,true,Fgt);

%%
% nBhv = 1;
A = randn(nBhv,sum(nF,'all')); 
if howManyColsPhi < 6
    A(:,howManyColsPhi+1:end) = 0; 
end
simulatedPsi = A * 1./vecnorm(A);
simulatedPsi(isnan(simulatedPsi))=0;

justOneGTStateExtra = groundTruthStatesExtra;
behaviorDataExtra = cellfun(@(x) simulatedPsi*x,justOneGTState,'UniformOutput',false);
%%
dFFExtra = statesToObservations(simulatedDmatrix, latentStatesX, noiseLevel);
save(sprintf('./TMLR/saveFMT_%d_seed09_TMLR_b%d_c%d_regen.mat',todaysDate,nBhv,howManyColsPhi))
%%
rng('shuffle')
currentState10 = rng;
%%
load(sprintf('./AISTATS/saveFMT_251124_seed10_AISTATS_b%d_c%d.mat',nBhv,howManyColsPhi))

[latentStatesXExtra, FgtNew, groundTruthStatesExtra] = generateContMultiSubNetwork_Range_Regenerate(...
               sizeT,[sizeD/2,sizeD/2],nF,50,contOpt,false,false,true,Fgt);

%%
% nBhv = 1;
A = randn(nBhv,sum(nF,'all')); 
if howManyColsPhi < 6
    A(:,howManyColsPhi+1:end) = 0; 
end
simulatedPsi = A * 1./vecnorm(A);
simulatedPsi(isnan(simulatedPsi))=0;

justOneGTStateExtra = groundTruthStatesExtra;
behaviorDataExtra = cellfun(@(x) simulatedPsi*x,justOneGTState,'UniformOutput',false);
%%
dFFExtra = statesToObservations(simulatedDmatrix, latentStatesX, noiseLevel);
save(sprintf('./TMLR/saveFMT_%d_seed10_TMLR_b%d_c%d_regen.mat',todaysDate,nBhv,howManyColsPhi))