function T = getPatternTable(Patterns)
% Takes pattern struct from a TheScript run and generates the corresponding
% table for that Pattern struct

% Label the first dimension of each struct
Patterns = nd.dimLabel(Patterns, 1, ["iPartition"]);

T = table();
for pattern = Patterns(:)'
    iPartition = pattern.iPartition;
    patternType = pattern.name;
    directionality = pattern.directionality;
    names = split(directionality, '-');
    source = names(1);
    target = names(2);
    nSource = size(pattern.X_source,1);
    nTarget = size(pattern.X_target,1);

    %Output properties
    rrDim = pattern.rankRegress.optDimReducedRankRegress;
    if isempty(rrDim)
        rankRegressDim = nan;
    end
    sourceFADim = pattern.factorAnalysis.qOpt_source;
    targetFADim = pattern.factorAnalysis.qOpt_target;
    if isempty(sourceFADim)
        sourceFADim = nan;
    end
    if isempty(targetFADim)
        targetFADim = nan;
    end
    maxDim = min(nSource,nTarget);
    percMax_rrDim = rrDim/maxDim;
    percMax_sourceFADim = sourceFADim/maxDim;
    percMax_targetFADim = targetFADim/maxDim;

    row = table(source, target, iPartition, patternType, nSource, nTarget, directionality, rrDim, sourceFADim, targetFADim, percMax_rrDim, percMax_sourceFADim, percMax_targetFADim);
    T = [T; row];
end
