function equalized = equalizePatternControl(nPatterns, cellOfWindows)
% equalized the number of windows between a pattern and its control
% cellOfWindows has to be formatted so that the first nPatterns are the
% patterns and the second half are the controls

equalized = cellOfWindows;

for pattern = 1:nPatterns
    curr_pattern = cell2mat(cellOfWindows(:,pattern));
    curr_control = cell2mat(cellOfWindows(:,pattern+nPatterns));
    [real_pattern, real_control] = ...
        windows.removeUntilEqual(curr_pattern, curr_control);
    equalized{pattern} = real_pattern;
    equalized{pattern+nPatterns} = real_control;
end

end

