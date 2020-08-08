function [H, times] = generateFromFilteredEEG(avgEEGStruct, brainArea, varargin)

% this function genegrates event windows based on averaged eeg files taking
% into consideration specified phase windows

% avgEEGStruct: eeg file containing filtered and averaged eeg
% brainArea: the area interested in (CA1)

% outputs:
% H: the event matrix (with specified rhythm patterns)

%% Parse optional arguments
ip = inputParser;
ip.addParameter('phaseWindow', []); % Window of possible phases to accept (if not empty)
ip.addParameter('downsample',  []); % How much (if not empty) to downsample data
ip.addParameter('patterns',    ["theta","delta","ripple"]); % Which patterns in avgEEGStruct to use
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

%% Iteratitvely build H
H = [];
times = [];
nPatterns = numel(opt.patterns);
for iPattern = 1:nPatterns
    patternH     = [];
  
    for index = indices'
        
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
        if ~isempty(opt.phaseWindow) % IF A PHASE WINDOW IS GIVEN BY THE USER, USE IT
            toNan = ~eventMatrix.inPhaseWindows(double(eegData(:,PHASE))/10000, ...
                opt.phaseWindow);
            eegData(toNan,:)   = nan;
        end
        
        patternH = [patternH;     eegData(:,AMP)];
        times    = [times; singleTime(:)];
    end
    % Build the data
    H     = [H patternH]; 
end

times = times';
