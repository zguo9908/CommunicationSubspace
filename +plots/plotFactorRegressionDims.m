function [outputArg1,outputArg2] = plotFactorRegressionDims...
                          (numDimsUsedForPrediction, cvLoss, qOpt varargin)
%PLOTFACTORREGRESSIONDIMS Summary of this function goes here
%   Detailed explanation goes here

ip = inputParser;
ip.addParameter('color', "black", @isstring);
ip.addParameter('averaged', true, @islogical);
ip.parse(varargin{:});
opt = ip.Results;

%% Cross-validate Factor Regression
% Plot Reduced Rank Regression cross-validation results
x = numDimsUsedForPrediction;
x(x > qOpt) = [];
y = 1-cvLoss(1,:);
e = cvLoss(2,:);

hold on
errorbar(x, y, e, 'o--', 'Color', COLOR(V2,:), ...
    'MarkerFaceColor', 'w', 'MarkerSize', 10)
hold off

end

