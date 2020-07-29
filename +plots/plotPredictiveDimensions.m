function plotPredictiveDimensions(numDimsUsedForPrediction, cvLoss, varargin)
% Plot Reduced Rank Regression cross-validation results

ip = inputParser;
ip.addParameter('color', "black",@isstring);
ip.parse(varargin{:});
opt = ip.Results;

x = numDimsUsedForPrediction;
y = 1-cvLoss(1,:);
e = cvLoss(2,:);

errorbarplot = errorbar(x, y, e);

errorbarplot.Color = opt.color;

xlabel('Number of predictive dimensions')
ylabel('Predictive performance')

end