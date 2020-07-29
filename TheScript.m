clear
%% Paths
if ~exist('animal','var')
    animal = "JS15";
end
addpath(genpath(pwd)) % all folders in utils added, including semedo code

% Whose computer are we running?
if ispc
    paths = "\\citadel.bio.brandeis.edu\sharespace-commsub\data\";
elseif ismac
    paths(1) = "/Volumes/sharespace-commsub/data";
    paths(2) = "~/Data/commsubspace";
end
% Set data paths for that computer
arrayfun(@(path) addpath(genpath(path)), paths);

%%
TABLE = load("Megatable.mat");

%% Script parameters
% -----------------------------------------------------------------
Option = struct();

Option.generateFromRipTimes = true; % Whether to replace H(:, ripple) as
% determined by spectral power with
% pattern determined by global ripple
% rimes.

Option.generateH = "fromFilteredEEG";
% this option field have the options 1)"fromSpectral", 2)"fromRipTimes",
% 3)"fromFilteredEEG"

Option.spikeBinSize = 0.1;          % 100 milliseconds
Option.timesPerTrial = 10;         % 10 times per trial
Option.winSize = 0.3;              % size of the window
Option.directionality = "2 direction";  %directionality (source-target)
Option.equalWindowsAcrossPatterns = false;    % whether all three patterns have the same #windows
Option.singleControl = true;                 % whether to use just one control column

usingSingleprediction = true;
%% Shortcut/alias variables to improve readability
THETA = 1;
DELTA = 2;
RIPPLE = 3;
HPC = 1;
PFC = 2;

%% Mung/Clean the data

load(animal + "spectralBehavior.mat");
spectrogram   = efizz.S1;
frequencyAxis = efizz.f;
times = efizz.t;
frequenciesPerPattern = [6 14; 0.5 4; 150 200];
[nPatterns,~] = size(frequenciesPerPattern);

if Option.singleControl == true
    numResult = nPatterns+1;
else
    numResult = nPatterns*2;
end

% find network pattern events
if contains(Option.generateH, "fromSpectral")
    H = eventMatrix.generateFromSpectra(times, spectrogram, frequencyAxis,...
                                             frequenciesPerPattern);
elseif contains(Option.generateH, "fromFilteredEEG")
    load(animal + "avgeeg.mat");
    [H, times] = eventMatrix.generateFromFilteredEEG(avgEEG, 'CA1')
end
if contains(Option.generateH,"fromRipTimes")
    load(animal + "globalripple01.mat");
    [ripplecolumntimes,ripplecolumn] = eventMatrix.generateFromRipples(globalripple, 500, ...
        'amplitude_at_riptime', true);
    samplepoints = interp1(ripplecolumntimes, ripplecolumn, times);
    H(:,RIPPLE) = samplepoints';
end

%%
%%%%%%%%%%%%%%%% WINDOW SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% windows of network patterns
cellOfWindows = windows.make(   times, 0.9, H(:,THETA:DELTA), Option.winSize);
cellOfWindows(3) = windows.make(times, 1,   H(:,RIPPLE),      Option.winSize, ...
    'threshold', 'raw');

% equalize number of windows across patterns based on input argument
if Option.equalWindowsAcrossPatterns == true
    cellOfWindows = windows.equalizeWindowsAcrossPatterns(cellOfWindows, nPatterns);
end

%% CONTROL SECTION
% add control patterns
Hc = control.generatePatternShuffle(H(:,1:3), times, cellOfWindows);

% add windows of control patterns
Hc_cellOfWindows = windows.make(times,  0.9, Hc(:,THETA:DELTA), Option.winSize);
Hc_cellOfWindows(3) = windows.make(times, 1, Hc(:,RIPPLE), Option.winSize, ...
    'threshold', 'raw');

% clean up control windows: remove each control pattern's window's overlap
for pattern = 1:nPatterns
    curr = windows.removeOverlapsBetweenPattern(...
        cell2mat(cellOfWindows(:,pattern)), ...
        cell2mat(Hc_cellOfWindows(:,pattern)));
    Hc_cellOfWindows{pattern} = curr;
end

% Merge into one
H(:,nPatterns+1:nPatterns*2) = Hc;
cellOfWindows(nPatterns+1:nPatterns*2) = Hc_cellOfWindows;

% Equalize trials/windows for each pair of patttern-controlPattern
cellOfWindows = control.equalizePatternControl(nPatterns, cellOfWindows);

%% Getting spikes
[timeAxis, times_spiking, spikeCountMatrix, spikeRateMatrix, areaPerNeuron] = ...
    spikes.getSpikeTrain(animal + "spikes01.mat", Option.spikeBinSize);

% making pattern matrices
[spikeSampleMatrix, spikeSampleTensor] = ...
    trialSpikes.generate(spikeCountMatrix, timeAxis, cellOfWindows, Option.timesPerTrial, numResult);

%% Separate spikesSampleMatrix/Tensor by area that neurons are in PFC and neurons that in HPC
X_pfc = trialSpikes.separateSpikes(spikeSampleMatrix, areaPerNeuron, "PFC");
X_hpc = trialSpikes.separateSpikes(spikeSampleMatrix, areaPerNeuron, "CA1");

%% Separate firing pattern into source and target
[nPFCneurons,~] = size(X_pfc{1});
[nHPCneurons,~] = size(X_hpc{1});

[X_source, X_target, nSource] = trialSpikes.splitSourceTarget(2, numResult, X_hpc, X_pfc, 'specifiedSource', "HPC");

%% Results place to store outputs
clear Patterns
Patterns.rankRegress = struct("B", [],"B_",[],...
    "optDimReducedRankRegress",0,"cvl_rankRegress",[],...
    "cvLoss_rankRegress",[],"singlesource_B",[],...
    "singlesource_optDim",[]);
Patterns.factorAnalysis = struct("qOpt_source",[],"Z_source",[],...
    "U_source",[],"Q_source",[],"qOpt_target",[],"Z_target",[],...
    "U_target",[],"Q_target",[]);

patternNames = ["theta","delta","ripple",...
    "theta-control","delta-control","ripple-control"];

for i = 1:numResult-1
    Patterns = [Patterns Patterns(1)];
end

if Option.directionality == "2 direction" % hpc-hpc and hpc-pfc
    Patterns = [Patterns; Patterns];
    for i = 1:numResult % "expected ONE OUTPUT from dot indexing"
        Patterns(HPC,i).directionality = "hpc-hpc";
        Patterns(HPC,i).X_source = X_source{i};
        Patterns(HPC,i).X_target = X_target{1,i};
        Patterns(HPC,i).name = patternNames(i);
        
        
        Patterns(PFC,i).directionality = "hpc-pfc";
        Patterns(PFC,i).X_source = X_source{i};
        Patterns(PFC,i).X_target = X_target{2,i};
        Patterns(PFC,i).name = patternNames(i);
    end
end

%% Regression cross-validation: FIND DIMENSION
numDimsUsedForPrediction = 1:nPFCneurons;

for i = 1:numResult
    
    B_singleprediction = cell(1,nSource);
    dim_singleprediction = cell(1,nSource);
    
    % Number of cross validation folds.
    cvNumFolds = 10;
    cvOptions = statset('crossval');
    regressMethod = @ReducedRankRegress;
    cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
        (regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
        numDimsUsedForPrediction, 'LossMeasure', 'NSE');
    
    curr_source = (Patterns(1,i).X_source)';
    
    for j = [HPC, PFC]
        curr_target = (X_target{j,i})';
        [Patterns(j,i).rankRegress.cvl, Patterns(j,i).rankRegress.cvLoss,...
            Patterns(j,i).rankRegress.optDimReducedRankRegress,Patterns(j,i).rankRegress.B,...
            Patterns(j,i).rankRegress.B_,Patterns(j,i).rankRegress.V] ...
            = rankRegressRoutine(cvFun, cvNumFolds, cvOptions,curr_target, curr_source, numDimsUsedForPrediction);
        try
            % Single neuron prediction
            for k = 1:nSource
                curr_singlesource = curr_source(:,j);
                [~,~, dim_singleprediction{k}, B_singleprediction{k},~,~] = ...
                    rankRegressRoutine(cvFun, cvNumFolds, cvOptions,curr_target, curr_singlesource,...
                    numDimsUsedForPrediction);
            end
            Patterns(j,i).rankRegress.singlesource_B = B_singleprediction;
            Patterns(j,i).rankRegress.singlesource_optDim = dim_singleprediction;
        catch
            usingSingleprediction = false;
        end
    end
    
end


%% Factor Analysis

doFactorAnalysis = true;
reset(gpuDevice()) %if the GPU has stuff in memory, clear the device, so we have a blank slate!
if doFactorAnalysis
    
    % How many dimensions did we maximally have for Bdims
    RR=[Patterns.rankRegress];
    maxBdim = max([RR.optDimReducedRankRegress]);
    
    %  Find factor analytic dimension
    cvNumFolds = 10;
    cvOptions = statset('crossval');
    cvOptions.UseParallel = true;
    fa_gpuState = true;
    if cvOptions.UseParallel && fa_gpuState
        parpool('local', 2); %  4   workers with gpu  crashes
    end
    
    
    for i = progress(1:1, 'Title', 'FA over patterns')
        
        
        curr_source = X_source{i};
        
        
        curr_pfc = X_pfc{i};
        
        q_pfc = 1:max(nPFCneurons,maxBdim);
        q_hpc = 1:max(nHPCneurons,maxBdim);
        
        tic % start timing how long this takes
        cvLoss_pfc= CrossValFa(curr_pfc, q_pfc, cvNumFolds, cvOptions, ...
            'gpu', fa_gpuState, 'speedup', false);
        disp(cvLoss_pfc);
        toc % print how long it took
        
        qOpt_pfc = FactorAnalysisModelSelect(cvLoss, q_pfc);
        
        tic;
        cvLoss_hpc= CrossValFa(curr_hpc, q_hpc, cvNumFolds, cvOptions, ...
            'gpu', fa_gpuState, 'speedup', false);
        qOpt_hpc = FactorAnalysisModelSelect(cvLoss, q_hpc);
        toc;
        
        Pattern(1,i).factorAnalysis.qOpt_hpc = qOpt_hpc;
        Pattern(2,i).factorAnalysis.qOpt_pfc = qOpt_pfc;
    end
    
    
    
    %     % -------------------------------------
    %     % For each area, extract latent factors
    %     % -------------------------------------
    %
    %     [Z_pfc, U_pfc, Q_pfc] = ExtractFaLatents(X_pfc', qOpt_pfc);
    %     [Z_hpc, U_hpc, Q_hpc] = ExtractFaLatents(X_hpc', qOpt_hpc);
end


%% Save results
optionFieldsNames =fieldnames(Option)';
Optiontable = table(optionFieldsNames);
Optiontable = table(struct2array(Option));
Patterntable = table(patternNames);
TABLE = table(Optiontable, Patterntable);

tableheader = table(optionFieldsNames, patternNames);
tablerow = table(struct2array(Option), Patterns(:)');
initialTable = [tableheader; tablerow];


TABLE = [TABLE.TABLE; tablerow];
save Megatable.mat TABLE


save('\\citadel.bio.brandeis.edu\sharespace-commsub/checkpoint_try.mat','-v7.3')

%% Higher priority TODO
% 2. Save table with all Option fields + number of hpc neurons useed to
% compute, number of pfc neurons used to compute, number of dimension found
% per pattern, cvLoss matrix (if possible), factor analysis outputs (we
% still lack these)
% 3. Implement factor analysis outputs
% 4. Implement acquiring Figure 2 distributions from semedo paper
% 5. Generalize window equalization more, at moment only option is to
% equalize numbers between pairs of windows

%% Lower priority TODO
% 2. Try Coherence Hs, SeqNMF Hs (I can supply), track spatial Hs

