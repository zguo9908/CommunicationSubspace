%function nonOverlappingpattern = removeOverlapWithinPattern(pattern)
    % pattern is a cell of start&stop times
    % returns a cell with start&stop times non-overlapping
    
    [a,nPatterns] = size(pattern);
    nonOverlappingpattern = cell(1,nPatterns);
    
    for i = 1:nPatterns
        
        patternWindows = pattern{i};
        [nWindows, ~] = size(patternWindows);
        result = patternWindows(1,:);
        
        tracker = 2; % a potential bug, now window 1 and 2 would ALWAYS be in the final list,
        % though they might be overlapping.
        while tracker < nWindows
            
            currWindow = patternWindows(tracker,:);
            nextWindow = patternWindows(tracker+1,:);
            
            if windows.detectOverlap(currWindow(1),currWindow(2), ...
                    nextWindow(1),nextWindow(2))
                
                patternWindows(tracker+1,:) = []; 
                
            end
            tracker = tracker+1;
            [nWindows, ~] = size(patternWindows);
            result = [result;currWindow];
             
        end
        
         nonOverlappingpattern{i} = result;
          
    end
   % nonOverlappingpattern = patternWindows;
    
%end