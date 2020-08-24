function [times, H, Hnanlocs, Hvals, valThreshold, original] = generateFromRipples(rippleData, varargin)
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
ip.addParameter('rippleBand',[]);
ip.addParameter('rippleBandTime',[]);
ip.addParameter('globalrippleWindowUnits','std');
ip.addParameter('samprate', 1500);
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
conditions = [isempty(opt.rippleBand), ...
              isempty(opt.rippleBandTime)];
if all(conditions)
    times = first_time_of_the_day:1/opt.samprate:last_time_of_the_day;
    H = interp1(opt.rippleBandTime,opt.rippleBand,times);
    Hvals = H; % this would create NaNs in Hvals
elseif any(conditions)
    error('Must provide both  rippleBand and original_time_axs');
else
    times = opt.rippleBandTime;
    H     = nan(1,length(times));
    switch opt.globalrippleWindowUnits
        case 'std'
            % Need to trasform rippleBand
            opt.rippleBand = (opt.rippleBand - mean(opt.rippleBand))./std(opt.rippleBand);
        case 'amp'
        otherwise
            error("Unhandled ripple unit");
    end
    Hvals = opt.rippleBand;
    original = true;
end


%% Create one list of ripples, all epochs
for iEpoch = 1:numel(ripples)
    if ~isempty(ripples{iEpoch})
        curr = ripples{iEpoch};
        for iWindows = 1:length(curr)
            ripple_windows = [ripple_windows;...
                              curr(iWindows,:)];
        end
    end
end

valThreshold = min(ripple_windows(:,3));

%%  Mark each time in the axis
for j = 1:size(ripple_windows,1)
    % first qualified window encountered
    filter = times >= ripple_windows(j,1) & times < ripple_windows(j,2);
    if opt.amplitude_at_riptime == false
        H(filter) = 1;
    else
        H(filter) = ripple_windows(j,3);
    end
    continue;
end

Hnanlocs = isnan(H);
Hnanlocs = double(Hnanlocs)';
Hnanlocs(Hnanlocs == 1) = nan;
Hnanlocs(Hnanlocs == 0) = 1;
Hvals(Hnanlocs==1) = H(Hnanlocs==1);

H = Hvals .* Hnanlocs;
