test = Patterns(1,1,1)
q = 0:40;

cvNumFolds_test = 10;

cvOptions = statset('crossval');
% cvOptions.UseParallel = true;
tic
cvLoss_test= CrossValFa(test.X_source', q, cvNumFolds_test, cvOptions);


% CrossValFa returns the cumulative shared variance explained. To compute
% the optimal Factor Analysis dimensionality, call
% FactorAnalysisModelSelect:
qOpt_test = FactorAnalysisModelSelect(cvLoss_test, q);


toc
%%
fa_gpuState = true;
if fa_gpuState
    reset(gpuDevice()) %if the GPU has stuff in memory, clear the device, so we have a blank slate!
end
if cvOptions.UseParallel && fa_gpuState
    parpool('local', 2); %  4   workers with gpu  crashes
end
    
tic
[Z, U, Q] = ExtractFaLatents(test.X_source', qOpt_test);
toc

%%
tic   
regressMethod = @FactorRegress;

p = size(test.X_source', 2);
qOpt = FactorAnalysisModelSelect( ...
	CrossValFa(test.X_source', q, cvNumFolds_test, cvOptions), ...
	q);
toc

tic
numDimsUsedForPrediction = 1:78;

cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
	(regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
	numDimsUsedForPrediction, ...
	'LossMeasure', 'NSE', 'qOpt', qOpt);
toc
%%
tic
cvl = crossval(cvFun, test.X_target', test.X_source', ...
	  'KFold', cvNumFolds_test, ...
	'Options', cvOptions);
toc

cvLoss = ...
	[ mean(cvl); std(cvl)/sqrt(cvNumFolds_test) ];
tic
optDimFactorRegress = ModelSelect...
	(cvLoss, numDimsUsedForPrediction);
toc