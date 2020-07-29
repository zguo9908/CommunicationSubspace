 function [cellOfWindows] = makeWindows(times, threshold, H, winSize, varargin)
% H: m*n matrix - m is number powers taken (aligned to times, same
%    length as times) and n is the number of rhythms we are interested in
%
% times: time points where the powers of different rhythms were taken
%
% threshold: percentile for which the powers are chosen
%
% winSize: length of window per oscillation
ip = inputParser;
ip.addParameter('thresholdMethod','quantile');
ip.parse(varargin{:})
opt = ip.Results;

%% 1 calculate the quantiles
[nTime,nPatterns] = size(H);
cutoff = zeros(1,nPatterns);
for i = 1:nPatterns
    % uncomment these for a histogram and the cutoff of the quantile
    %          figure(i); clf
    %          hcount = histogram(H(:,i)); % RY: What am I doing here? What's returned by histogram()
    switch opt.thresholdMethod
        case 'quantile'
        cutoff(i) = quantile(H(:,i),threshold);
        case 'raw'
        cutoff(i) = threshold;
        otherwise
            error('Invalid method')
    end

    %             hold on
    %           lineObject=line([cutoff(i),cutoff(i)],[0 max(hcount.Values)]);
    %     lineObject.LineStyle = ':'; % Make line dotted
    %     lineObject.LineWidth = 2;  % Thicken the line
    %     lineObject.Color = 'black'; % Color it black
    %     xlim([0 max(H(:,1))]);
    %     M = area([cutoff(i),max(H(:,1))],[max(hcount.Values) max(hcount.Values)],"FaceAlpha",0.5);
    %     M.FaceColor = [0.85,0.33,0.10];
end

%% 2 Apply quantiles to the powers

%marking the H_matrix with 1 if the power is above quantile and 0 if
%not
for i = 1:length(cutoff)
    for j = 1:nTime
        if H(j,i)>=cutoff(i)
            H(j,i)=1;
        else
            H(j,i)=0;
        end
    end
end

%% 3 map to the time points
cellOfWindows = cell(1,length(cutoff));
% the time points we want, 1st column starting and 2nd ending
for i = 1:length(cutoff)
    %result{i} = zeros(1,2);
    switches = diff(H(:,i));
    switches = (switches == 1);
    temp = find(switches);
    temp = (switches == 1);
    R = times(temp)';
    
    if isscalar(winSize)
        R(:,2) = R(:,1)+winSize;
    elseif numel(winSize) == 2
        R = [R+winSize(1), R+winSize(2)];
    else
        error("winSize either a single number (total window time from the trigger) or two numbers (time in front of and behind the trigger)")
    end
    
    cellOfWindows{i} = R;
end
end
