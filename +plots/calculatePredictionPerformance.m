function [performance, mean_performance, nan_indices] = ...
                   calculatePredicationPerformance(X_source, X_target, B)
% calculate the predictive performance of source firing to target firing

[~,nTarget] = size(X_target);
try
[loss, yhat] = RegressPredict(X_target, X_source, B);
catch
    keyboard
end
performance = [];

for j = 1:nTarget
    unaccounted_variance =sum((yhat(:,j)-X_target(:,j)).^2);
    total_variance = sum((mean(X_target(:,j))-X_target(:,j)).^2);
    
  
    temp_r_square = 1-unaccounted_variance/total_variance;
%     
%     if isnan(temp_r_square)
%         temp_r_square = -10000;
%     end
%     
    performance = [performance temp_r_square];
end

nan_indices = find(isnan(performance));

indices = intersect(find(~isinf(performance)), find(~isnan(performance)));

mean_performance = mean(performance(indices));
end

