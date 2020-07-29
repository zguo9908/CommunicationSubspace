function [times, ripples_are_happening] = generateFromRipples(rippleData, samprate, varargin)
% Input
% -----
% rippleData : cell of structs
%   branched cell of structs loaded from *ripples.mat
%
% samperate : scalar
%   the sampling rate you desire for your output
%
% (optional keyword)
%
% amplitude_at_riptime

%
% Output
% ------
% times : 1 x T
%     Literally, first_time_of_the_day : 1/samprate : last_time_of_the_day
% ripples_are_happening : 1 x T
%     At each time a record of whether ripples are happening


ip = inputParser;
ip.addParameter('amplitude_at_riptime', false, @islogical);
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
ripples_are_happening = zeros(1,length(times));
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
for i = 1:length(times)
    ripples_that_are_in_the_past = ripple_windows(:,2) < times(i);
    ripple_windows(ripples_that_are_in_the_past,:) = []; % RY: Delete ripples with ends less than the current time ... speeds up code, lessens search time
    for j = 1:size(ripple_windows,1)
        if times(i)>=ripple_windows(j,1) && times(i)<ripple_windows(j,2)
            if opt.amplitude_at_riptime == false 
                ripples_are_happening(i) = 1;
            else
                ripples_are_happening(i) = ripple_windows(j,3);
            end
            break;
        end
    end
end

end