% function controlWindows = generateControlWindows(cellOfWindows, times, winSize)
% -------------------------------------------------------------------------%
% input
% cellOfWindows: the cell of matices of start & stop times for each pattern
% times: the time axis
% winSize: size of window of high activity
%
% output 
% controlWindows: the cell of matrices of start & stop times 
%--------------------------------------------------------------------------%

%% deselect periods of high activity from times based on cellOfWindows
[a, nPatterns] = size(cellOfWindows);
possibleWindows = cell(1,nPatterns);
for i = 1:nPatterns
    curr_times = times;
    curr_pattern = cellOfWindows{i};
    [nWindows, ~] = size(curr_pattern);
    for j = 1:nWindows
        window_time = times>=curr_pattern(j,1)&times<curr_pattern(j,2);
        curr_times(window_time) = 0;
    end
    possibleWindows{i} = curr_times(curr_times~=0);
end

%% pick control windows by start times that are possible
controlWindows = cell(1,nPatterns);
for i = 1:nPatterns
    curr_possible = possibleWindows{i};
    curr_pattern = cellOfWindows{i};
    nWindows = length(curr_possible);
    curr_windows = [];
    for j = 1:nWindows
        curr_windows = [curr_windows; curr_possible(j) curr_possible(j)+winSize];
    end
    controlWindows{i} = curr_windows;
end

%% removeOverlapBetweenWindows per pattern and its matching control
nonOverlappingpattern = cell(1,nPatterns);

for i = 1:nPatterns
    curr_controlPattern = controlWindows{i};
    curr_Pattern = cellOfWindows{i};
    
    %test
%    curr_controlPattern = curr_controlPattern(1:10000,:);
%     curr_Pattern = curr_Pattern(1:200,:);
    
%not really removing
    nonOverlappingpattern{i} = removeOverlapsBetweenPattern(curr_Pattern, curr_controlPattern);
end

%% remove patterns for so that the control and the actual patterns have equal amount of windows

    %removeWindowsUntilEqual
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%     [nPossibletime, nPossiblePattern] = size(possible_times);
%     whole_control_column = zeros(nTimes,nPatterns); 
%     if nPossiblePattern == 3
%            
%         for i = 1:nPatterns
%             curr_power = H(:,i);
%             for j = 1:nTimes % similar comment here for line 40
%                 shuffled_sequence = randperm(nTimes);
%                 k = 2;
%                 while ~ismember(shuflled_sequence(1),possible_times{i})
%                     shuffle_sequence(1) = shuffled_sequence(k);
%                     k = k+1;
%                 end
%                 shuffled = curr_power(randperm(nTimes));
%                 
%                 whole_control_column(j,i) = curr(1);
%             end
%             toAdd = toAdd';
%             controlH = [controlH toAdd];
%         end
%         
%     else   
%         for iPattern = 1:nPatterns
%             %curr_power = H(:,randi(3));
%             curr_power = H(:,iPattern); % you just take the whole column. This is going to be permuted.
%             whole_control_column(:,iPattern) = curr_power(randperm(nTimes));            
%         end
%         controlH = [controlH whole_control_column];
%         
%     end
% end