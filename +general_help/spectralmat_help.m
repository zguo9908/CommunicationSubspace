% spectral.mat
%
% How can we obtain C, S1, S2?

% C : matrix of time x frequency 
%    PFC-HPC Coherenece in all frequencies over all times

% S1 : matrix of time x frequency 
%    HPC spectral power for all time and frequencies

% S2 : matrix of time x frequency 
%    PFC spectral power for all frequencies over all times

% t : time axis of the three above

% f : frequency axis of the three above
%
%

%% Example: Plotting small snippet of time in each

load('spectral.mat')

% First grab a subset of times from each
small_subset_of_time = 1000:2000;
t = t(small_subset_of_time);
C = C(small_subset_of_time, :);
S1 = S1(small_subset_of_time, :);
S2 = S2(small_subset_of_time, :);

% Grab a subset of frequencies, only 0-80hz let's says

% why 109 frequencies when f>0 and f<80
small_frequency_subset = f > 0 & f < 80;
f = f(small_frequency_subset);
C = C(:, small_frequency_subset);
S1 = S1(:, small_frequency_subset);
S2 = S2(:, small_frequency_subset);


% S2 and S1 probably need a good ol fashioned zscore to see much of
% anything
zS1 = zscore(S1, 0, 1); % Zscore S1 without sample size correction (0) and over time-axis (1)
zS2 = zscore(S2, 0, 1);

figure;
subplot(2,2,1)
imagesc(t, f, C'); set(gca,'ydir','normal'); % Transpose and flip each... same as doing imagesc(t,f,C'); set(gca,'ydir','normal');
xlabel("Time (Seconds)"); ylabel("Frequency (Hertz)")
title('PFC HPC Coherence')
colorbar;
subplot(2,2,2)
imagesc(t, f, zS2'); set(gca,'ydir','normal'); 
xlabel("Time (Seconds)"); ylabel("Frequency (Hertz)")
title('PFC Spectra')
colorbar;
subplot(2,2,3)
imagesc(t, f, zS1'); set(gca,'ydir','normal'); 
xlabel("Time (Seconds)"); ylabel("Frequency (Hertz)")
title('HPC Spectra')
colorbar;


% Or log10 scaling
% why would this be helpful, start from 1?
logS1 = log10(S1+10); % Zscore S1 without sample size correction (0) and over time-axis (1)
logS2 = log10(S2+10);

figure;
subplot(2,2,1)
imagesc(t, f, C'); set(gca,'ydir','normal'); % Transpose and flip each... same as doing imagesc(t,f,C'); set(gca,'ydir','normal');
xlabel("Time (Seconds)"); ylabel("Frequency (Hertz)")
title('PFC HPC Coherence')
subplot(2,2,2)
imagesc(t, f, logS2'); set(gca,'ydir','normal'); 
xlabel("Time (Seconds)"); ylabel("Frequency (Hertz)")
title('PFC Spectra')
subplot(2,2,3)
imagesc(t, f, logS1'); set(gca,'ydir','normal'); 
xlabel("Time (Seconds)"); ylabel("Frequency (Hertz)")
title('HPC Spectra')
