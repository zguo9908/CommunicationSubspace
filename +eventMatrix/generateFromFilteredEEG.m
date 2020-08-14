function [H, Hvals, Hnanlocs, times] = generateFromFilteredEEG(avgEEGStruct, task, brainArea, varargin)

% this function genegrates event windows based on averaged eeg files taking
% into consideration specified phase windows

% avgEEGStruct: eeg file containing filtered and averaged eeg
% task: the taskfile information so that sleep sessions can be ignored
% brainArea: the area interested in (CA1)

% outputs:
% H: the event matrix (with specified rhythm patterns)
% Hvals: version of H where instead of nanning out the rhythm not
% happening, their actual values are used
% Hnanlocs: where the nans are in the final H matrix
% times: the time axis

%% Parse optional arguments
ip = inputParser;
ip.addParameter('phaseWindow', [0-pi/4 0+pi/4]); % Window of possible phases to accept (if not empty)
ip.addParameter('downsample',  []); % How much (if not empty) to downsample data
ip.addParameter('patterns',    ["theta","delta","ripple"]); % Which patterns in avgEEGStruct to use
ip.addParameter('ignoreSleep', true);
ip.parse(varargin{:});
opt = ip.Results;

PHASE = 2; AMP = 3;

% If user passes a branched cell array, instead of an ndimension struct, converet it!
if iscell(avgEEGStruct)
    indices = ndBranch.indicesMatrixForm(avgEEGStruct); %get the indices of all the cells within the avg eeg files, each with pfc and hpc info
    avgEEGStruct = ndBranch.toNd(avgEEGStruct); %convert into structs with relevant information
else
    indices = nd.indicesMatrixForm(avgEEGStruct);
end


%% get epochs with running sessions

sleepSessions = [];

if opt.ignoreSleep 
    task = task{1};
    
    for i = 1:numel(task)
        if task{i}.type == "sleep"
            sleepSessions = [sleepSessions, i];
        end
    end
% else, there will be no sleepSessions, and all the structs will be
% examined.
end 

%% Iteratitvely build H
Hvals = [];
Hnanlocs = [];
nPatterns = numel(opt.patterns);

for iPattern = 1:nPatterns
    patternH     = [];
    patternToNan = [];
    times = []; % Do not move this  line: If placed above for iPattern = 1:nPatterns, creates a bug
   
    for index = indices'
        if ismember(index(2), sleepSessions)
            continue
        end
        
        I = num2cell(index);
        singleStruct = avgEEGStruct( I{:} );
        
        % If this data not from the requested brain area, skip
        if ~isequal( singleStruct.area , brainArea )
            continue
        end
        
        % Get times of this single day-ep-brainarea
        singleTime = geteegtimes(singleStruct);
        % Data for ths single  day-ep-brainarea
        eegData     = singleStruct.(opt.patterns(iPattern)).data;
        eegData = double(eegData);
        
        % DOWNSAMPLE?
        if ~isempty(opt.downsample) % IF A DOWNSAMPLE IS GIVEN BY THE USER, USE IT
             eegData     = downsample(eegData,opt.downsample);
             singleTime = downsample(singleTime,opt.downsample);
        end

        % SELECT PHASES?
        toNan = zeros(1,length(eegData));
        if ~isempty(opt.phaseWindow) % IF A PHASE WINDOW IS GIVEN BY THE USER, USE IT
            toNan = ~eventMatrix.inPhaseWindows(double(eegData(:,PHASE))/10000, ...
                opt.phaseWindow);
%             eegData(toNan,:)   = nan;
        end
        
        patternH = [patternH;     eegData(:,AMP)];
        patternToNan = [patternToNan;  toNan];
        times    = [times; singleTime(:)];
    end
    
    patternNanlocs = ones(1,length(patternToNan));
    patternNanlocs(patternToNan==1) = nan;
        
    % Build the data
    Hvals     = [Hvals patternH]; 
    Hnanlocs  = [Hnanlocs patternNanlocs'];
end

H =  Hvals .* Hnanlocs;

times = times';
