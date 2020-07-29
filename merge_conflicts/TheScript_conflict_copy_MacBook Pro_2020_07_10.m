
clear

%% Paths
if ~exist('animal','var')
    animal = "JS15";
end
addpath(genpath("data" + filesep + animal)); % genpath enumerates every folder in inside the folder it's given
addpath(genpath('utils')); % all folders in utils added, including semedo code
TABLE = readtable("Megatable.csv");

% -----------------------------------------------------------------
%% Script parameters
% -----------------------------------------------------------------
Option = struct();

Option.generateFromRipTimes = false; % Whether to replace H(:, ripple) as 
                                    % determined by spectral power with
                                    % pattern determined by global ripple
                                    % rimes.
                                    
Option.spikeBinSize = 0.1;          % 100 milliseconds
Option.timesPerTrial = 10;         % 10 times per trial
Option.winSize = 0.25;              % size of the window
Option.directionality = "hpc-pfc";  %directionality (source-target)
Option.compute.fa_gpuState = true;

%% Shortcut/alias variables to improve readability
THETA = 1;
DELTA = 2;
RIPPLE = 3;

%% Mung/Clean the data
%change here for different animals 
load(animal + "spectralBehavior.mat");
spectrogram   = efizz.S1;
frequencyAxis = efizz.f;
times = efizz.t;
frequenciesPerPattern = [6 14; 0.5 4; 150 200];
[nPatterns,~] = size(frequenciesPerPattern);

% find network pattern events
H = eventMatrix.generateFromSpectra(times, spectrogram, frequencyAxis,...
    frequenciesPerPattern); % (A)
if Option.generateFromRipTimes == true
    load(animal + "globalripple01.mat");

    [ripplecolumntimes,ripplecolumn] = eventMatrix.generateFromRipples(globalripple, 500, ...
        'amplitude_at_riptime', true); % RY look at your output here, compare to the output of the function file
    %interval = length(ripplecolumn)/length(H(:,RIPPLE));
    %samplepoints = ripplecolumn(1:interval:end);
    samplepoints = interp1(ripplecolumntimes, ripplecolumn, times);
    H(:,RIPPLE) = samplepoints';
end
% % windows of network patterns
[cellOfWindows] = windows.make(times, 0.9, H, Option.winSize); %(B)


%%%%%%%%%%%%%%%% CONTROL SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add control patterns
Hc = control.generatePatternShuffle(H, times, cellOfWindows); %(A)
% add windows of control patterns
Hc_cellOfWindows = windows.make(times, 0.9, Hc, Option.winSize);

% clean up control windows: remove each control pattern's window's overlap
for pattern = 1:nPatterns
    curr = windows.removeOverlapsBetweenPattern(...
        cell2mat(cellOfWindows(:,pattern)), ...
        cell2mat(Hc_cellOfWindows(:,pattern)));
    Hc_cellOfWindows{pattern} = curr;
end

% Merge into one
H = [H, Hc];
cellOfWindows = [cellOfWindows, Hc_cellOfWindows];

% Equalize trials/windows for each pair of patttern-controlPattern
for pattern = 1:nPatterns
    curr_pattern = cell2mat(cellOfWindows(:,pattern));
    curr_control = cell2mat(cellOfWindows(:,pattern+nPatterns));
    [real_pattern, real_control] = ...
        windows.removeUntilEqual(curr_pattern,curr_control);
    cellOfWindows{pattern} = real_pattern;
    cellOfWindows{pattern+nPatterns} = real_control;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting spikes
[timeAxis, times_spiking, spikeCountMatrix, spikeRateMatrix, areaPerNeuron] = ...
    spikes.getSpikeTrain(animal + "spikes01.mat", Option.spikeBinSize);

% making pattern matrices
[spikeSampleMatrix, spikeSampleTensor] = ...
    trialSpikes.generate(spikeCountMatrix, timeAxis, cellOfWindows, Option.timesPerTrial);

for pattern = 1:nPatterns*2
    curr = spikeSampleMatrix{pattern};
    curr(isnan(curr)) = 0;
    spikeSampleMatrix{pattern} = curr;
end

%% Separate spikesSampleMatrix/Tensor by area that neurons are in PFC and neurons that in HPC
X_pfc = trialSpikes.separateSpikes(spikeSampleMatrix, areaPerNeuron, "PFC");
X_hpc = trialSpikes.separateSpikes(spikeSampleMatrix, areaPerNeuron, "CA1");

%% Results struct
clear Patterns
Patterns = struct("name","","X_pfc",[],"X_hpc",[],"B_",[],...
    "optDimReducedRankRegress",0,"cvl_rankRegress",[],...
    "cvLoss_rankRegress",[],"qOpt_pfc",[],"Z_pfc",[],"U_pfc",[],"Q_pfc",[],...
    "qOpt_hpc",[],"Z_hpc",[],"U_hpc",[],"Q_hpc",[]);

patternNames = ["theta","delta","ripple",...
    "theta-control","delta-control","ripple-control"];
for i = 1:nPatterns*2-1
    Patterns = [Patterns Patterns(1)];
end

labels = ["\theta", "\delta", "\rho"];
labels = [labels; labels + "-control"]';
for i = 1:nPatterns*2
    Patterns(i).X_pfc = X_pfc{i};
    Patterns(i).X_hpc = X_hpc{i};
    Patterns(i).label = labels(i);
    Patterns(i).name = patternNames(i);
end

[nPFCneurons,~] = size(X_pfc{1});
[nHPCneurons,~] = size(X_hpc{1});

%% Regression cross-validation: FIND DIMENSION ++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% % Cross-validate Reduced Rank Regression
for i = 1:nPatterns*2
    X = X_hpc{i};
    Y_V2 = X_pfc{i};

    % Vector containing the interaction dimensionalities to use when fitting
    % RRR. 0 predictive dimensions results in using the mean for prediction.
    numDimsUsedForPrediction = 1:nPFCneurons;

    % Number of cross validation folds.
    cvNumFolds = 10;

    % Initialize default options for cross-validation.
    cvOptions = statset('crossval');

    regressMethod = @ReducedRankRegress;

    cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
        (regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
        numDimsUsedForPrediction, 'LossMeasure', 'NSE');

    % Cross-validation routine.
    cvl = crossval(cvFun, Y_V2', X', ...
          'KFold', cvNumFolds, ...
        'Options', cvOptions);

    % Stores cross-validation results: mean loss and standard error of the
    % mean across folds.
    cvLoss = [ mean(cvl); std(cvl)/sqrt(cvNumFolds) ];

    % To compute the optimal dimensionality for the regression model, call
    % ModelSelect:
    optDimReducedRankRegress = ModelSelect...
        (cvLoss, numDimsUsedForPrediction);

    Patterns(i).optDimReducedRankRegress = optDimReducedRankRegress;
    Patterns(i).cvl = cvl;
    Patterns(i).cvLoss = cvLoss;
end

%%  Reduced rank regression w/ dimension
% % -------------------------------------
for i = 1:nPatterns*2
    curr_hpc = X_hpc{i};
    curr_pfc = X_pfc{i};
    optimalDim = Patterns(i).optDimReducedRankRegress;
    [B, B_, V] = ReducedRankRegress(curr_pfc', curr_hpc', optimalDim);
    Patterns(i).B_ = B_;
    Patterns(i).B = B;
    Patterns(i).V = V;
    %[B, B_, V] = ReducedRankRegress(X_pfc',X_hpc', 1:nPFCneurons, 'RIDGEINIT', true);
end


%% Factor Analysis

%  Find factor analytic dimension


cvNumFolds = 10;
% k
cvOptions = statset('crossval');
cvOptions.UseParallel = false;
% 
for i = progress(1:nPatterns*2)
    curr_hpc = X_hpc{i};
    curr_pfc = X_pfc{i};
    
    if Option.compute.fa_gpuState
        [curr_hpc, curr_pfc] = deal(gpuArray(curr_hpc), ...
                                    gpuArray(curr_pfc));
    end
    
    [nPFCneurons, nHPCneurons] = deal(size(curr_hpc,1),size(curr_pfc,1));  
    
    q_pfc = 1:nPFCneurons;
    q_hpc = 0:nHPCneurons;
    
    tic % start timing how long this takes
    cvLoss_pfc= CrossValFa(curr_pfc, q_pfc, cvNumFolds, cvOptions, ...
        'gpu', Option.compute.fa_gpuState);
    disp(cvLoss_pfc);
    toc % print how long it took

    qOpt_pfc = FactorAnalysisModelSelect(cvLoss, q_pfc);
        
%     cvLoss_hpc= CrossValFa(curr_hpc, q_hpc, cvNumFolds, cvOptions);
% 
%     qOpt_hpc = FactorAnalysisModelSelect(cvLoss, q_hpc);
end



% -------------------------------------
% For each area, extract latent factors
% -------------------------------------

% [Z_pfc, U_pfc, Q_pfc] = ExtractFaLatents(X_pfc', qOpt_pfc);
% [Z_hpc, U_hpc, Q_hpc] = ExtractFaLatents(X_hpc', qOpt_hpc);


%% Save results

% Save what  you've computed thus  far  to a TABLE to with parameters of interests and to 
% checkpoint file.
dimensions = extractfield(Patterns,'optDimReducedRankRegress');
newrow = {Option.generateFromRipTimes, Option.spikeBinSize, Option.timesPerTrial, ...
    Option.winSize, char(Option.directionality),nPFCneurons,nHPCneurons,mat2str(dimensions)};
TABLE = [TABLE;newrow];
writetable(TABLE,'Megatable.csv');

save('../checkpoint.mat','-v7.3') 

%% Higher priority TODO
% 1. Visualize windows on raw data
% 2. Save table with all Option fields + number of hpc neurons useed to
% compute, number of pfc neurons used to compute, number of dimension found
% per pattern, cvLoss matrix (if possible), factor analysis outputs (we
% still lack these)
% 3. Implement factor analysis outputs
% 4. Implement acquiring Figure 2 distributions from semedo paper

% 5. Generalize window equalization more, at moment only option is to
% equalize numbers between pairs of windows
% 6

%% Lower priority TODO
% 1. Option.directionality = "hpc-hpc", "hpc-pfc" etc
% 2. Try Coherence Hs, SeqNMF Hs (I can supply), track spatial Hs

%% ----------------------Notes---------------------------------------------
% 1. Maybe ridge or scale is not the best
% 2. Maybe windows need to be improved
%   a. Window size (milliseconds per window or number of time points per window)
%        1. Total width of window
%        2. Total number of points
%   b. Threshold that triggers : quantile threshold to trigger a pattern
%       1. Too low? Too high?
%       2. Noise sometimes in lfp can create windows which have random cell
%       activity
% 3. Maybe we need more neurons : use ZT2 dataset (justin's 64 tetrode
% animal on citadel or thunderdome).

