clear
%% Paths
addpath(genpath(pwd)) % all folders in utils added, including semedo code


if ~exist('animal','var')
    animal = "ZT2";
end
addpath(genpath(pwd)) % all folders in utils added, including semedo code

% Whose computer are we running?
if ispc
    paths = " C:\Users\BrainMaker\commsubspace\SingleDayExpt";
elseif ismac
    paths(1) = "/Volumes/sharespace-commsub/data";
    paths(2) = "~/Data/commsubspace";
end
% Set data paths for that computer
arrayfun(@(path) addpath(genpath(path)), paths);
% % ---------- Add paths for specific users ------------------------------
% addpath(genpath(datadef(user))); % this line throws an error for pc
%
% addpath(...
%     fullfile(codedefine,"Shared"));
% addpath(genpath(...
%     fullfile(codedefine, "Shared", "utils")...
%     ));
% % -----------------------------------------------------------------------
%% Script parameters
% -----------------------------------------------------------------
Option = struct();

Option.generateFromRipTimes = true; % Whether to replace H(:, ripple) as
% determined by spectral power with
% pattern determined by global ripple
% rimes.
Option.animal    = animal;
Option.generateH = "fromSpectral_fromRipTimes";
%Option.generateH = "fromSpectra_fromRipTimes";
% This option field have these possibilities:
%  1)"fromSpectral",
%  2)"fromFilteredEEG",
%  3)"fromFilteredEEG_fromRipTimes"
%  4)"fromSpectra_fromRipTimes", where fromRipTimes causes only the ripple
%  pattern to  derive from global ripples

Option.spikeBinSize  = 0.1;          % 100 milliseconds
Option.timesPerTrial = 10;         % 10 times per trial
Option.winSize       = 0.3;              % size of the window
Option.sourceArea    = "CA1";
Option.equalWindowsAcrossPatterns = false;    % whether all three patterns have the same #windows
Option.singleControl = true;                 % whether to use just one control column
Option.numPartition = 10;                    % ways to split source and target
usingSingleprediction = true;
%% Shortcut/alias variables to improve readability
THETA = 1;
DELTA = 2;
RIPPLE = 3;
HPC = 1;
PFC = 2;

%% Mung/Clean the data
frequenciesPerPattern = [6 14; 0.5 4; 150 200];
[nPatterns,~] = size(frequenciesPerPattern);

if Option.singleControl == true
    numResult = nPatterns+1;
else
    numResult = nPatterns*2;
end

% find network pattern events
if contains(Option.generateH, "fromSpectral")
    load(Option.animal + "spectralBehavior.mat");
    if Option.sourceArea == "CA1"
        spectrogram   = efizz.S1;
    else
        spectrogram   = efizz.S2;
    end
    
    frequencyAxis = efizz.f;
    times = efizz.t;
    H = eventMatrix.generateFromSpectra(times, spectrogram, frequencyAxis,...
        frequenciesPerPattern);
elseif contains(Option.generateH, "fromFilteredEEG")
    load(Option.animal + "avgeeg.mat");
    [H, times] = eventMatrix.generateFromFilteredEEG(avgeeg, ...
        Option.sourceArea, "patterns",["theta","delta","ripple"],"downsample",100);
end

%making the ripple column of H zero when H is generated from Filtered EEG?
if contains(Option.generateH,"fromRipTimes")
    load(Option.animal + "globalripple01.mat");
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
[cellOfWindows, warnedEmptyControls] = control.equalizePatternControl(nPatterns, cellOfWindows);

% pick control pattern that actually contains controls, would break if all
% three are empty...
if Option.singleControl
    if warnedEmptyControls
        for iPossibleControl = nPatterns+1:nPatterns*2
            if ~isempty(cellOfWindows{iPossibleControl})
                cellOfWindows{nPatterns+1} = cellOfWindows{iPossibleControl};
                disp("here")
            end
        end
    end
end

%%%%%%%%%%%%%%%% SPIKE SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Getting spikes
[timeAxis, times_spiking, spikeCountMatrix, spikeRateMatrix, areaPerNeuron] = ...
    spikes.getSpikeTrain(Option.animal, Option.spikeBinSize);

% making pattern matrices
[spikeSampleMatrix, spikeSampleTensor] = ...
    trialSpikes.generate(spikeCountMatrix, timeAxis, cellOfWindows, Option.timesPerTrial, numResult);

%% Separate spikesSampleMatrix/Tensor by area that neurons are in PFC and neurons that in HPC
X_pfc = trialSpikes.separateSpikes(spikeSampleMatrix, areaPerNeuron, "PFC");
X_hpc = trialSpikes.separateSpikes(spikeSampleMatrix, areaPerNeuron, "CA1");

%% Separate firing pattern into source and target
[nPFCneurons,~] = size(X_pfc{1});
[nHPCneurons,~] = size(X_hpc{1});

X_source = cell(Option.numPartition, 1, numResult);
X_target = cell(Option.numPartition, 2, numResult);
% 
% [X_source, X_target, nSource] = trialSpikes.splitSourceTarget...
%     (2, numResult, X_hpc, X_pfc, 'specifiedSource', "CA1");

%%%%%%%%%%%%%%%% PRE-RESULT SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results place to store outputs
clear Patterns
Patterns = struct("X_source",[], "X_target",[]);
Patterns.rankRegress = struct(...
    "B", [], ...
    "B_", [], ...
    "optDimReducedRankRegParess", 0, ...
    "cvl_rankRegress", [], ...
    "cvLoss_rankRegress", [], ...
    "singlesource_B", [], ...
    "singlesource_optDim",[]);
Patterns.factorAnalysis = struct(...
    "qOpt_source",[],...
    "Z_source",[],...
    "U_source",[],...
    "Q_source",[],...
    "qOpt_target",[],...
    "Z_target",[],...
    "U_target",[],...
    "Q_target",[]);

patternNames = ["theta","delta","ripple",...
    "theta-control","delta-control","ripple-control"];


Patterns = repmat(Patterns, [Option.numPartition,2,numResult]);
%
% for i = 1:numResult-1
%     Patterns = [Patterns Patterns(1)];
% end
%
%
% Patterns = [Patterns; Patterns];
for iPartition = 1:Option.numPartition
    [X_source(iPartition,:),X_target(iPartition,:,:), nSource] = trialSpikes.splitSourceTarget...
        (2, numResult, X_hpc, X_pfc, 'specifiedSource', "CA1");
    
    for i = 1:numResult % "expected ONE OUTPUT from dot indexing"
        Patterns(iPartition,HPC,i).directionality = "hpc-hpc";
        Patterns(iPartition,HPC,i).X_source = X_source{iPartition,:,i};
        Patterns(iPartition,HPC,i).X_target = X_target{iPartition,1,i};
        Patterns(iPartition,HPC,i).name = patternNames(i);
        
        
        Patterns(iPartition,PFC,i).directionality = "hpc-pfc";
        Patterns(iPartition,PFC,i).X_source = X_source{iPartition,:,i};
        Patterns(iPartition,PFC,i).X_target = X_target{iPartition,2,i};
        Patterns(iPartition,PFC,i).name = patternNames(i);
    end
end
%%
%%%%%%%%%%%%%%%% RANK-REGRESS SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numDimsUsedForPrediction = 1:nPFCneurons;
for p = 1:Option.numPartition
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
    
    curr_source = (Patterns(p,1,i).X_source)';
    
    for j = [HPC, PFC]
        
        curr_target = (X_target{p,j,i})';
        [   Patterns(p,j,i).rankRegress.cvl, ...
            Patterns(p,j,i).rankRegress.cvLoss,...
            Patterns(p,j,i).rankRegress.optDimReducedRankRegress,...
            Patterns(p,j,i).rankRegress.B,...
            Patterns(p,j,i).rankRegress.B_,...
            Patterns(p,j,i).rankRegress.V] ...
            = rankRegressRoutine(cvFun, cvNumFolds, ...
            cvOptions, curr_target, curr_source, ...
            numDimsUsedForPrediction);
        try
            % Single neuron prediction
            for k = 1:nSource
                curr_singlesource = curr_source(:,j);
                if clean.zeroFiring(curr_singlesource)
                    continue;
                end
                [~,~, ...
                    dim_singleprediction{k}, ...
                    B_singleprediction{k},~,~] = ...
                    rankRegressRoutine(cvFun, cvNumFolds, ...
                    cvOptions,curr_target, ...
                    curr_singlesource,...
                    numDimsUsedForPrediction);
            end
            Patterns(p,j,i).rankRegress.singlesource_B = B_singleprediction;
            
            Patterns(p,j,i).rankRegress.singlesource_optDim = ...
                dim_singleprediction;
        catch
            usingSingleprediction = false;
        end
    end
    
end
end












%%
%%%%%%%%%%%%%%%% FACTOR ANALYSIS SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doFactorAnalysis = false;
if doFactorAnalysis
    
    % How many dimensions did we maximally have for Bdims
    RR=[Patterns.rankRegress];
    maxBdim = max([RR.optDimReducedRankRegress]);
    minBdim = min([RR.optDimReducedRankRegress]);
    
    %  Find factor analytic dimension
    cvNumFolds = 10;
    cvOptions = statset('crossval');
    cvOptions.UseParallel = true;
    fa_gpuState = true;
    if fa_gpuState
        reset(gpuDevice()) %if the GPU has stuff in memory, clear the device, so we have a blank slate!
    end
    if cvOptions.UseParallel && fa_gpuState
        parpool('local', 2); %  4   workers with gpu  crashes
    end
    
    
    for i = progress(1:nPatterns, 'Title', 'FA over patterns')
        
        
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
%%

%%%%%%%%%%%%%%%% CREATE TABLE AND SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare table!
TABLE = load("Megatable.mat");
Optiontable  = array2table(struct2array(Option));
Patterntable = array2table( Patterns(:)');

% Identifying information about this options set and date of run
hash = DataHash(Option);
hash = hash(1:7); % Take the first 7 letters of the hash
hash = string(hash);
datestr   = cell(1,1);
datastr{1} = date();
Filetable = table(hash, datestr);

% If it's a previous seen option, replace the row, otherwise append
tablerow = table(Optiontable, Patterntable, Filetable);
% Prevous options: Replace row
%%
if any(contains(TABLE.Filetable.hash, hash))
    TABLE(contains(TABLE.Filetable.hash, hash), :) = tablerow;
    disp("already computed before, rehashing to the same location");
    % New options:    Append row
else
    TABLE    = [TABLE.TABLE; tablerow];
    disp("new results stored!")
end

%%
% Save
save Megatable.mat TABLE
save(fullfile(datadefine, "hash", hash), ...
    "Patterns", "Option", '-v7.3')

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

