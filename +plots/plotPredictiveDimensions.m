function plotPredictiveDimensions(numDimsUsedForPrediction, cvLoss, varargin)
% Plot Reduced Rank Regression cross-validation results

ip = inputParser;
ip.addParameter('color', "black",@isstring);
ip.parse(varargin{:});
opt = ip.Results;

x = 1:numDimsUsedForPrediction;
y = 1-cvLoss(1,:);
e = cvLoss(2,:);
try
errorbarplot = errorbar(x, y, e);
catch
    keyboard
end
errorbarplot.Color = opt.color;

xlabel('Number of predictive dimensions')
ylabel('Predictive performance')

end