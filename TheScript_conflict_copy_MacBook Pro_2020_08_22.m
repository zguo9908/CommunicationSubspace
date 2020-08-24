% clear
% %% Paths
% addpath(genpath(pwd)) % all folders in utils added, including semedo code
% 
% 
% if ~exist('animal','var')
%     animal = "JS17";
% end
% addpath(genpath(pwd)) % all folders in utils added, including semedo code
% 
% % Whose computer are we running?
% if ispc
%     paths = " C:\Users\BrainMaker\commsubspace\SingleDayExpt";
% elseif ismac
%     paths(1) = "/Volumes/sharespace-commsub/data";
%     paths(2) = "~/Data/commsubspace";
% end
% % Set data paths for that computer
% arrayfun(@(path) addpath(genpath(path)), paths);
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
% % -----------------------------------------------------------------
% Option = struct();

Option.generateFromRipTimes = true; % Whether to replace H(:, ripple) as
% determined by spectral power with
% pattern determined by global ripple
% rimes.
% Option.animal    = animal;
Option.generateH = "fromFilteredEEG "+" fromRipTimes";
% Option.generateH = "fromSpectral "+" fromRipTimes";
% Option.generateH = "fromCoherence "+" fromRipTimes";
% This option field have these possibilities:
%  1)"fromSpectral",
%  2)"fromFilteredEEG",
%  pattern to derive from global ripples

Option.samplingRate  = nan ;              % For spikes.getSpikeTrain, nan if not given
Option.spikeBinSize  = 0.1;               % 100 milliseconds
Option.timesPerTrial = 10;                % 10 times per trial
Option.winSize       = {[-0.2, 0.2]};     % size of the window
Option.sourceArea    = "CA1";
Option.equalWindowsAcrossPatterns = true;    % whether all three patterns have the same #windows
Option.singleControl = false;                 % whether to use just one control column
Option.numPartition = 10;                    % ways to split source and target
usingSingleprediction = true;

winSize = Option.winSize{1};
%% Shortcut/alias variables to improve readability
THETA = 1;
DELTA = 2;
RIPPLE = 3;
if Option.sourceArea == "CA1"
    HPC = 1;
    PFC = 2;
else
    PFC = 1;
    HPC = 2;
end
patternNames = ["theta","delta","ripple"];

%% Mung/Clean the data
frequenciesPerPattern = [6 10; 0.5 4; 150 200];
[nPatterns,~] = size(frequenciesPerPattern);

if Option.singleControl == true
    numResult = nPatterns+1;
else
    numResult = nPatterns*2;
end

% Find network pattern events
if contains(Option.generateH, "fromSpectral")
    load(Option.animal + "spectralBehavior.mat");
    if Option.sourceArea == "CA1"
        spectrogram   = efizz.S1;
    else
        spectrogram   = efizz.S2;
    end
    
    frequencyAxis = efizz.f;
    Htimes = efizz.t;
    [H, Hvals, Hnanlocs, Htimes] = eventMatrix.generateFromSpectra(Htimes, spectrogram, frequencyAxis,...
        frequenciesPerPattern);
elseif contains(Option.generateH, "fromFilteredEEG")
    load(Option.animal + "avgeeg.mat");
    [~, sleepSessions] = getRunningSessions(Option.animal);
    [H, Hvals, Hnanlocs, Htimes] = eventMatrix.generateFromFilteredEEG(avgeeg, ...
        Option.sourceArea, "patterns",patternNames(1:3),"downsample",10, "sleepSessions", sleepSessions);
elseif contains(Option.generateH, "fromCoherence")
   load(Option.animal + "spectralBehavior.mat");
   spectrogram = efizz.C;
    frequencyAxis = efizz.f;
    Htimes = efizz.t;
    [H, Hvals, Hnanlocs, Htimes] = eventMatrix.generateFromSpectra(Htimes, spectrogram, frequencyAxis,...
        frequenciesPerPattern);
end
%%

if contains(Option.generateH,"fromRipTimes")
    load(Option.animal + "globalripple01.mat");
   
    [~, H(:,RIPPLE), Hnanlocs(:,RIPPLE), Hvals(:,RIPPLE), minRippleThreshold, original] = ...
        eventMatrix.generateFromRipples(globalripple, ...
                                       'amplitude_at_riptime', true,...
    
end
%%
%%%%%%%%%%%%%%%% WINDOW SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% windows of network patterns
cellOfWindows = windows.make(   Htimes, 0.9, H(:,THETA:DELTA), winSize);
cellOfWindows(3) = windows.make(Htimes, 1,   H(:,RIPPLE),      winSize, ...
    'threshold', 'raw');

disp(newline);
disp('Initiale window  creation:')
disp('--------------------------')
cellfun(@(x,y)  fprintf("%d windows for %s\n", size(x,1),  y), ...
    cellOfWindows, cellstr(patternNames(1:3)));

% equalize number of windows across patterns based on input argument
if Option.equalWindowsAcrossPatterns == true
    cellOfWindows = windows.equalizeWindowsAcrossPatterns(cellOfWindows, nPatterns);
end

disp(newline);
disp('Equalized windows        :')
disp('--------------------------')
cellfun(@(x,y)  fprintf("%d windows for %s\n", size(x,1),  y), ...
    cellOfWindows, cellstr(patternNames(1:3)));

%windows.printWindowOverlap(cellOfWindows, patternNames);

%% CONTROL SECTION
% add control patterns
Hc = control.generatePatternShuffle(H(:,1:3), Htimes, cellOfWindows);

% add windows of control patterns
Hc_cellOfWindows = windows.make(Htimes,  0.9, Hc(:,THETA:DELTA), winSize);
Hc_cellOfWindows(3) = windows.make(Htimes, 1, Hc(:,RIPPLE), winSize, ...
    'threshold', 'raw');

% clean up control windows: remove each control pattern's window's overlap
for pattern = 1:nPatterns
    curr = windows.removeOverlapsBetweenPattern(...
     cell2mat(cellOfWindows(:,pattern)), cell2mat(Hc_cellOfWindows(:,pattern)));
    Hc_cellOfWindows{pattern} = curr;
end


% % Merge into one
% H(:,nPatterns+1:nPatterns*2) = Hc;
cellOfWindows(nPatterns+1:nPatterns*2) = Hc_cellOfWindows;

% Equalize trials/windows for each pair of patttern-controlPattern
[cellOfWindows, warnedEmptyControls] =...
    control.equalizePatternControl(nPatterns, cellOfWindows);

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

if ~any(contains(patternNames,"control"))
    patternNames = [patternNames; patternNames+"-control"]';
    patternNames = patternNames(:)';
end

disp(newline);
disp('Equalized control windows :')
disp('--------------------------')
cellfun(@(x,y)  fprintf("%d windows for %s\n", size(x,1),  y), ...
    cellOfWindows, cellstr(patternNames));

%%%%%%%%%%%%%%%% SPIKE SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Getting spikes
[timeAxis, times_spiking, spikeCountMatrix, spikeRateMatrix, areaPerNeuron] = ...
    spikes.getSpikeTrain(Option.animal, Option.spikeBinSize,  Option.samplingRate);

% making pattern matrices
[spikeSampleMatrix, spikeSampleTensor] = ...
    trialSpikes.generate(spikeCountMatrix, timeAxis, cellOfWindows, ...
    Option.timesPerTrial, numResult);

%% Separate spikesSampleMatrix/Tensor by area that neurons are in PFC and neurons that in HPC
X_pfc = trialSpikes.separateSpikes(spikeSampleMatrix, areaPerNeuron, "PFC");
X_hpc = trialSpikes.separateSpikes(spikeSampleMatrix, areaPerNeuron, "CA1");

%% Separate firing pattern into source and target
[nPFCneurons,~] = size(X_pfc{1});
[nHPCneurons,~] = size(X_hpc{1});

X_source = cell(Option.numPartition, 1, numResult);
X_target = cell(Option.numPartition, 2, numResult);

%%%%%%%%%%%%%%%% PRE-RESULT SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results place to store outputs
clear Patterns
Patterns = struct("X_source",[], "X_target",[]);
Patterns.rankRegress = struct(...
    "B", [], ...
    "B_", [], ...
    "optDimReducedRankRegress", 0, ...
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

if numResult == nPatterns+1
    patternNames = ["theta","delta","ripple","control"];
end


Patterns = repmat(Patterns, [Option.numPartition,2,numResult]);

for iPartition = 1:Option.numPartition
    [X_source(iPartition,:),X_target(iPartition,:,:), nSource, nTarget] ...
        = trialSpikes.splitSourceTarget...
        (2, numResult, X_hpc, X_pfc, 'specifiedSource', "CA1");
    
    for i = 1:numResult % "expected ONE OUTPUT from dot indexing"
        
        if Option.sourceArea == "CA1"
            Patterns(iPartition,HPC,i).directionality = "hpc-hpc";
            Patterns(iPartition,PFC,i).directionality = "hpc-pfc";
        else
            Patterns(iPartition,HPC,i).directionality = "pfc-hpc";
            Patterns(iPartition,PFC,i).directionality = "pfc-pfc";
        end
        Patterns(iPartition,HPC,i).X_source = X_source{iPartition,:,i};
        Patterns(iPartition,HPC,i).X_target = X_target{iPartition,1,i};
        Patterns(iPartition,HPC,i).name = patternNames(i);
        
        Patterns(iPartition,PFC,i).X_source = X_source{iPartition,:,i};
        Patterns(iPartition,PFC,i).X_target = X_target{iPartition,2,i};
        Patterns(iPartition,PFC,i).name = patternNames(i);
    end
end
%%
%%%%%%%%%%%%%%%% RANK-REGRESS SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numDimsUsedForPrediction = 1:nTarget;
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
    q = 1: nTarget-1;
    cvNumFolds = 10;
    cvOptions = statset('crossval');
    for p = 1:Option.numPartition
        for i = 1:numResult
            for j = [HPC, PFC]
                disp("processing the "+p+" partition and the "+i+" pattern"+j)
                currSource = Patterns(p,j,i).X_source';
                cvLoss = CrossValFa(currSource, q, cvNumFolds, cvOptions);
                qOpt =  FactorAnalysisModelSelect(cvLoss, q);
                Patterns(p,j,i).factorAnalysis.cvLoss = cvLoss;
                Patterns(p,j,i).factorAnalysis.qOpt = qOpt
                [Patterns(p,j,i).factorAnalysis.Z, ...
                 Patterns(p,j,i).factorAnalysis.U,...
                 Patterns(p,j,i).factorAnalysis.Q] = ...
                 ExtractFaLatents(X_source', qOpt);
            end
        end
    end
end
%%

%%%%%%%%%%%%%%%% CREATE TABLE AND SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare table!
path2 = "C:\Users\BrainMaker\MATLAB Drive\Shared";
cd(path2)
load("Megatable.mat");
winSizes = cell2mat(Option.winSize);
if numel(Option.winSize{1} == 2)
    % otherwise, could not put a cell in a table, and cell2mat creates
    % varying length of vectors which will be counted as
    Option.winSize = "pre: "+winSizes(1)+ " post: "+winSizes(2);
else
    Option.winSize =  " post: "+winSizes(1);
end

Optiontable  = struct2table(Option);
Patterntable = query.getPatternTable(Patterns);

% Identifying information about this options set and date of run
hash = DataHash(Option);
hash = hash(1:7); % Take the first 7 letters of the hash
hash = string(hash);
timestamp = string(datestr());
tablerow = [Optiontable, table(timestamp, hash)]; % Combine option columnns with hash and date
tablerow = [Patterntable, repmat(tablerow, height(Patterntable), 1)]; % combine those with all rows of the pattern table

%% Check and Hash
if any(contains(TABLE.hash, hash))
    TABLE(contains(TABLE.hash, hash), :) = []; % Delete any rows that contain the current hash
    disp("already computed before, rehashing to the same location");
    % New options:    Append row
else
    TABLE    = [TABLE; tablerow];
    disp("new results stored!")
end

%% Save
save ("Megatable", "TABLE",'-v7.3');
save(fullfile(datadefine, "hash", hash), ...
    "Patterns", "Option", "H", "cellOfWindows", '-v7.3')
%%
% 2. Save table with all factor analysis outputs (we
% still lack these)
% 3. Implement factor analysis outputs
% 4. Implement acquiring Figure 2 distributions from semedo paper

