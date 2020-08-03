function inPhase = inPhaseWindows(phases, phaseWindows)
% this function examines if a particular phase belong to a set of phase
% windows, which are passed in as a single phase window 1*2 matrix.
    
inPhase =        phases <= phaseWindows(2) & ...
                 phases > phaseWindows(1) & ~isnan(phases);
