function [times, H, Hnanlocs, Hvals, original] = generateFromRipples_new...
    (rippleData, ripple_to_replace, original_time, varargin)

ip = inputParser;

ip.addParameter('samplingMethod', "original", @isstring);
% the new ripple column can be sampled either by a given sampling rate or
% by the rate of the ripple column passed in
ip.addParameter('samprate', []);
% if chosen to sample at a certain rate ("customized"), need to pass it in

ip.parse(varargin{:});
opt = ip.Results;

session = 1;
ripples = rippleData{session};

ripple_windows = [];

%% Create one list of ripples, all epochs
for iEpoch = 1:numel(ripples)
    if ~isempty(ripples{iEpoch})
        curr = ripples{iEpoch};
        for iWindows = 1:length(curr)
            ripple_windows = [ripple_windows;curr(iWindows,:)];
        end
    end
end

%% filling in the ripple column
Hvals = ripple_to_replace;
if opt.samplingMethod == "original"
    original = true;
    Hvaluelocs = [];
    for i = 1:length(original_time)
        
        ripples_that_are_in_the_past = ripple_windows(:,2) < original_time(i);
        % stop searching when the end time surpasses the time point interested in
        ripple_windows(ripples_that_are_in_the_past,:) = [];
        for j = 1:size(ripple_windows,1)
            % first qualified window encountered
            if  original_time(i)>=ripple_windows(j,1) && ...
                    original_time(i) < ripple_windows(j,2)
                Hvaluelocs = [Hvaluelocs, i];               
                Hvals(i) = ripple_windows(j,3);
                break;
            end
        end
    end
    
    H = Hvals;
    H(~Hvaluelocs) = nan;
    Hnanlocs = isnan(H);
    times = original_time;
else
    % how to fill in the Hvals for times that are sampled? 
    [first_time_of_the_day, last_time_of_the_day] = getTimerange(ripples);
%     first_time_of_the_day = intmax;
%     last_time_of_the_day = -1;
%     % find the inital and end time of the day
%     for iEpoch = 1:numel(ripples)
%         if ~isempty(ripples{iEpoch})
%             curr = ripples{iEpoch};
%             curr_first = curr(1,1);
%             curr_last = curr(end,2);
%             if curr_first<first_time_of_the_day
%                 first_time_of_the_day = curr_first;
%             end
%             if curr_last>last_time_of_the_day
%                 last_time_of_the_day = curr_last;
%             end
%         end
%     end
    
    times = first_time_of_the_day:1/samprate:last_time_of_the_day;
    
    temp = interp1(original_time, ripple_to_replace, times);
    
    for i = 1:length(times)
        ripples_that_are_in_the_past = ripple_windows(:,2) < times(i);
        % stop searching when the end time surpasses the time point interested in
        ripple_windows(ripples_that_are_in_the_past,:) = [];
        for j = 1:size(ripple_windows,1)
            % first qualified window encountered
            if times(i)>=ripple_windows(j,1) && times(i) < ripple_windows(j,2)
                    H(i) = ripple_windows(j,3);
                break;
            end
        end
    end
    Hvaluelocs = isnan(H); % where all the ripple events don't occur
    Hvals(~Hvaluelocs) = H(~Hvaluelocs);
end

