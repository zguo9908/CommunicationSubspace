function plotPredictiveDimensions(numDimsUsedForPrediction, cvLoss, varargin)
% Plot Reduced Rank Regression cross-validation results

ip = inputParser;
ip.addParameter('color', "black", @isstring);
ip.addParameter('averaged', true, @islogical);
ip.parse(varargin{:});
opt = ip.Results;

x = 1:numDimsUsedForPrediction;

if ~opt.averaged
    y = 1-cvLoss(1,:);
    e = cvLoss(2,:);
    
else
    numPartitions = numel(cvLoss);
    toAverage = [];
    error = [];
    for i = 1:numPartitions
        curr_cvLoss = cvLoss{i};
        toAverage = [toAverage; curr_cvLoss(1,:)];
        error = [error; curr_cvLoss(2,:)];
    end
    y = 1-mean(toAverage,1);
    e = mean(error,1);

end


    errorbarplot = errorbar(x, y, e);

errorbarplot.Color = opt.color;

xlabel('Number of predictive dimensions')
ylabel('Predictive performance')

end