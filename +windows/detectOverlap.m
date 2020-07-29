function isOverlapping = detectOverlap(start1, end1, start2, end2)
    % detect if two windows overlap

    if start1>end1 || start2>end2
        error("start time must be smaller than end time for both patterns.");
    end
    
    if start1<=end2 && start2<=end1
        isOverlapping = true;
    else
        isOverlapping = false;
    end
    
end