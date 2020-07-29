function [H, times] = generateFromFilteredEEG(avgEEGStruct, brainArea, varargin)
%GENERATEFROMFILTEREDEEG (in devo, not tested, possibly flawed)
% this function genegrates event windows based on averaged eeg files taking
% into consideration specified phase windows

% avgEEGStruct:

%% Parse optional arguments
ip = inputParser;
ip.addParameter('phaseWindow', []); % Window of possible phases to accept (if not empty)
ip.addParameter('downsample',  []); % How much (if not empty) to downsample data
ip.addParameter('patterns',    ["theta","delta"]); % Which patterns in avgEEGStruct to use
ip.parse(varargin{:});
opt = ip.Results;

% If user passes a branched cell array, instead of an ndimension struct, converet it!
if iscell(avgEEGStruct)
    indices = ndBranch.indicesMatrixForm(avgEEGStruct); %get the indices of all the cells within the avg eeg files, each with pfc and hpc info
    avgEEGStruct = ndBranch.toNd(avgEEGStruct); %convert into structs with relevant information
else
    indices = nd.indicesMatrixForm(avgEEGStruct);
end

%% Iteratitvely build H
times = [];
H     = [];

for index = indices'
  count = count+1
    I = num2cell(index);
    singleStruct = avgEEGStruct( I{:} );

        % If this data not from the requested brain area, skip
        if ~isequal( singleStruct.area , brainArea )
            continue
        end
        
        % Get times of this single day-ep-brainarea
        singleTimes = geteegtimes(singleStruct);
        
        % Get the H of this single day-ep-brainarea
        % (This set of two steps might be wrong if singleH cell have structs)
        singleH    = arrayfun(@(field) singleStruct.(field), opt.patterns,...
            'UniformOutput', false);
%          singleH    = arrayfun(@(field) singleStruct.(field), "delta",...
%             'UniformOutput', false);
        singleH     = cat(2, singleH{1}.data); 
    
    % first downsample 
    if ~isempty(opt.downsample)
        singleH = downsample(singleH, opt.downsample);
        singleTimes = downsample(singleTimes, opt.downsample);
    end
%     singleH = downsample(singleH,100);
%         singleTimes = downsample(singleTimes,100);
    if ~isempty(opt.phaseWindow)
         todelete = [];
        for iPossibleWindow = 1:size(singleH,1)
            if (~eventMatrix.inPhaseWindows(double(singleH(iPossibleWindow,2))/10000,opt.phaseWindow))
                todelete = [todelete iPossibleWindow];
            end
        end
        singleH(todelete,:)   = []; % RY: nan or zero maybe instead of delete to cut down on discontinuity?.. especially when you plot this, it may look a little funky even though correct.
        singleTimes(todelete) = []; % and if nans/zeros, then comment this out.
        % nans probably over zeros ... nans do not contribute to the
        % quantile stats during the windows.make() call
    end

    % Build the data
    times = [times singleTimes];
    H     = [H;     singleH(:,1)];

end
