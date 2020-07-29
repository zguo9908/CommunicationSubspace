function [controlH] = generateControlColumn(H, times, exclusionTimesPerPattern)


[nTime, nPattern] = size(H);

%% Make a shuffled iPattern
% -------------------------
controlH = H;
for iPattern = 1:nPattern
    controlH(:,iPattern) = H(randperm(nTime),iPattern); 
    % shuffle the current H matrix
end

%% Exclude times
% --------------
% To exclude times, we just set event matrix entries to zero 
% and they will not be selected for windows downstream
START = 1; % Just a constant/alias to improve readability below
STOP  = 2;
for iPattern = 1:nPattern
    nWindows       = size(exclusionTimesPerPattern,1);
    exclusionTimes = exclusionTimesPerPattern{iPattern};
    for iWindow = 1:nWindows
        timesToZeroOut = times > exclusionTimes(iWindow, START ) & ...
            times <= exclusionTimes(iWindow, STOP);
        controlH(timesToZeroOut, iWindow) = 0;
    end
end

