function [cvl, cvLoss, optDim, B, B_, V] = rankRegressRoutine(cvFun, cvNumFolds, cvOptions,target, source, num)
% repetitive steps when running rank regress

cvl = crossval(cvFun, target, source, ...
    'KFold', cvNumFolds, ...
    'Options', cvOptions);
cvLoss = [ mean(cvl); std(cvl)/sqrt(cvNumFolds)];
optDim = ModelSelect...
    (cvLoss, num);

[B, B_, V] = ReducedRankRegress(target,source,optDim);


end

