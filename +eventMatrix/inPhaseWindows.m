function inPhase = inPhaseWindows(singlePhase, phaseWindows)
% this function examines if a particular phase belong to a set of phase
% windows, which are passed in as n*2 matrix.

inPhase = false;
nPossiblePhases = size(phaseWindows, 1);
for i = 1:nPossiblePhases
    
    if ~isnan(singlePhase)
        
        if singlePhase < phaseWindows(i,2) && singlePhase>phaseWindows(i,1)
            inPhase = true;
        end
    end
end

end

