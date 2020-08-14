function [times, H, Hnanlocs, Hvals] = generateFromRipples(rippleData, ...
    samprate, varargin)
% Input
% -----
% rippleData : cell of structs
%   branched cell of structs loaded from *ripples.mat
%
% samperate : scalar
%   the sampling rate you desire for your output
%
% (optional keyword) amplitude_at_riptime: if true, return the amplitute

%
% Output
% ------
% times : 1 x T
%     Literally, first_time_of_the_day : 1/samprate : last_time_of_the_day
% H: the matrix where non-rippling times are marked nans while the rippling
%    events are marked with their ripple power if chose for amplitute, or 1 if
%     not.
% Hnanlocs: where in the H column the nans are
% Hvals: all the ripple power of the H column


ip = inputParser;
ip.addParameter('amplitude_at_riptime', false, @islogical);
ip.addParameter('ripple_to_replace',[]);
ip.addParameter('original_time_axis',[]);
ip.parse(varargin{:});
opt = ip.Results;

session = 1;
ripples = rippleData{session};
first_time_of_the_day = intmax;
last_time_of_the_day = -1;
ripple_windows = [];

%% Find the initial and end times of the day
for iEpoch = 1:numel(ripples)
    if ~isempty(ripples{iEpoch})
        curr = ripples{iEpoch};
        curr_first = curr(1,1);
        curr_last = curr(end,2);
        if curr_first<first_time_of_the_day
            first_time_of_the_day = curr_first;
        end
        if curr_last>last_time_of_the_day
            last_time_of_the_day = curr_last;
        end
    end
end

%% Creat time axis and initialize output
times = first_time_of_the_day:1/samprate:last_time_of_the_day;
if ~isempty(opt.ripple_to_replace) && ~isempty(opt.original_time_axis)
    H = interp1(opt.original_time_axis,opt.ripple_to_replace,times);
    Hvals = H;
else
    H = nan(1,length(times));
    Hvals = zeros(1,length(times));
end


%% Create one list of ripples, all epochs
for iEpoch = 1:numel(ripples)
    if ~isempty(ripples{iEpoch})
        curr = ripples{iEpoch};
        for iWindows = 1:length(curr)
            ripple_windows = [ripple_windows;curr(iWindows,:)];
        end
    end
end

%%  Mark each time in the axis
for j = 1:size(ripple_windows,1)
    %ripples_that_are_in_the_past = ripple_windows(:,2) < times(i);
    % stop searching when the end time surpasses the time point interested in
    %ripple_windows(ripples_that_are_in_the_past,:) = [];
        % first qualified window encountered
    filter = times >= ripple_windows(j,1) & times < ripple_windows(j,2); % RY instead of for-looping times, we vectorize here for speed
    if any(filter)
        if opt.amplitude_at_riptime == false
            H(filter) = 1;
        else
            H(filter) = ripple_windows(j,3);
        end
        break;
    end
end

Hnanlocs = isnan(H); % where all the ripple events don't occur
Hvals(~Hnanlocs) = H(~Hnanlocs);
end
