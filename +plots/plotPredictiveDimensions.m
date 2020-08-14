function full_model = plotPredictiveDimensions(numDimsUsedForPrediction, cvLoss, varargin)

% Plot Reduced Rank Regression cross-validation results
% also returns the averaged full model
% full_model = 1-curr_cvLoss(1,end);
ip = inputParser;
ip.addParameter('color', "black", @isstring);
ip.addParameter('averaged', true, @islogical);
ip.parse(varargin{:});
opt = ip.Results;

x = 1:numDimsUsedForPrediction;

if ~opt.averaged
    y = 1-cvLoss(1,:);
    e = cvLoss(2,:);
    full_model = 1-cvLoss(1,end);
else
    numPartitions = numel(cvLoss);
    toAverage = [];
    error = [];
    all_full_model = [];
    for i = 1:numPartitions
        curr_cvLoss = cvLoss{i};
        toAverage = [toAverage; curr_cvLoss(1,:)];
        error = [error; curr_cvLoss(2,:)];
        all_full_model = [all_full_model, 1-curr_cvLoss(1,end)];
    end
    y = 1-mean(toAverage,1);
    e = mean(error,1);
    full_model = mean(all_full_model);
end


errorbarplot = errorbar(x, y, e);

errorbarplot.Color = opt.color;

xlabel('Number of predictive dimensions')
ylabel('Predictive performance')

end