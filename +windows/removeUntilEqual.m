function [pattern, control]=removeUntilEqual(curr_pattern, curr_control)
    
    % this function takes in two pattern windows and remove the one with
    % more windows until they are equal.
    
    [n1,~] = size(curr_pattern);
    [n2,~] = size(curr_control);
    
    if n1>n2
        toRemove = curr_pattern;
    else
        toRemove = curr_control;
    end
    removeSize = abs(n1-n2);
    
    indices = randperm(removeSize);
    toRemove(indices',:) = [];
 
    if n1>n2
        pattern = toRemove;
        control = curr_control;
    else
        control = toRemove;
        pattern = curr_pattern;
    end
end