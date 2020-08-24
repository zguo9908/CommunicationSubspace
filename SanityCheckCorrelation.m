

load(animal+"pos01.mat");
posdata = pos{1}{2}.data;
postime = unique(posdata(:,1));
posvel = posdata(:,end);
Hvals(:,RIPPLE) = H(:,RIPPLE);
Hvals(isnan(Hvals(:,RIPPLE)),RIPPLE) =0;
indices = BetweenTimes(Htimes, postime); % since H_times has the higher sampling rate, we want to translate the higher to the lower
indices = indices(~isnan(indices));

% posvel = interp1(postime, posvel, indices);
% Theta has positive speed correlation
velocity_correlations = corrcoef(...
[posvel, Hvals(indices,:)])
velocity_correlations(isnan(velocity_correlations))=0;

figure(145)
for i = 1:size(velocity_correlations,1); velocity_correlations(i,i)=nan; end
imagesc(velocity_correlations)
patternNames = ["velocity","theta","delta","ripple","theta-control","delta-control","ripple-control"];
set(gca,'XTick', [1:7],'XTickLabel',patternNames);
set(gca,'YTick',[1:7],"YTickLabel",patternNames);
cmocean('balance')
colorbar
set(gca,'clim',[-1,1])
title(plots.getOptionInfo(Option))