% To code your FUNCTION, try first writing it as a
% SCRIPT. Use the SCRIPT to try out changes to the code easily. Then after
% it your
% sure you're algorithm works. When it does, convert to a function!

% function [cellOfWindows] = makeWindows(times, Q, H, winSize) % Uncomment
% this line when your script works
% H: m*n matrix - m is number powers taken (aligned to times, same
%    length as times) and n is the number of rhythms we are interested in
% times: time points where the powers of different rhythms were taken
% Q: percentile for which the powers are chosen
% winSize: length of window per oscillation

%[newTimes, animalStartTimes] = makeSuperAnimal(times)

%% setting up inputs
clear
addpath('./Shared')
load('test-window-function')

% Input
Q = 0.9;
[times, startTimes] = makeSuperAnimal(t);
winSize = 100;
H = H';


%% 1 calculate the quantiles

num = size(H); % number of network patterns to find quantiles for
cutoff = zeros(1,num(2));
for i = 1:num(2)
    figure(i); clf
    hcount = histogram(H(:,i)); % RY: What am I doing here? What's returned by histogram()
    %ZG: the sizes of the bins (frequency)
    cutoff(i) = quantile(H(:,i),Q);
    hold on
    % Make vertical line
    % ------------------
    lineObject=line([cutoff(i),cutoff(i)],[0 max(hcount.Values)]); % RY: Why did I change this: 
    % (1) notice I'm giving two arguments to line(), not three
    % (2) was max(H(:,i)) the max of a given histogram? What is on the y-axis of a histogram? And where are values of H(:,i)? on the x? on the y?
    % Style the line object
    % ---------------------
    lineObject.LineStyle = ':'; % Make line dotted
    lineObject.LineWidth = 2;  % Thicken the line
    lineObject.Color = 'black'; % Color it black
    % Create a shaded area to the right of the vertical line
    % ------------------------------------------------------
    % How would you shade the area right of the dotted quantile y-line red?
    % (see area() command ... trying to teach figure skillz)
    xlim([0 max(H(:,1))]);
    M = area([cutoff(i),max(H(:,1))],[max(hcount.Values) max(hcount.Values)],"FaceAlpha",0.5);
    M.FaceColor = [0.85,0.33,0.10];
end
%% 2 Apply quantiles to the powers
for i = 1:length(cutoff)
    for j = 1:num(1)
        if H(j,i)>=cutoff(i)
            H(j,i)=1;
        else
            H(j,i)=0;
        end
    end
end

%% 3 map to the time points
result = cell(1,length(cutoff)); 
% the time points we want, 1st column starting and 2nd ending
for i = 1:length(cutoff)
    %result{i} = zeros(1,2);
    switches = diff(H(:,i)); % two issues: (see below)
    % (1) do you really want the absolute value here? what do the three possible values {0,1,-1} mean? in the command window below, try plotting H(:,1) with stairs() and then plot the switches as vertical lines were you think window starts are happening ...
          %ZG: I see, 0 means within a window, 1 means switching on, -1 is
          %off. However, as we've defined winSizes, we don't need the -1's
          %I guess. But in order to logically index, we could only have 0/1
          %
          %switches = switches(switches == 1); why this does not work to
          %filter out the 1's?
          % 
          %% RY: so that expression is just going to return all the 1s and no
          % 0s. To understand, look as these two expressions (and try to
          % understand why they're different)
          %
          % (Expression 1) switches = (switches == 1); % Correct ... 
          % returns 1s at each location of 1 and 0 for any other value, same length as switches
          %
          % (Expression 2) switches = switches(switches == 1); % Incorrect ... 
          % returns a sequences of pure 1s with length smaller than the original switches
          %
          % Expression 2 is the same thing as:
          % switches = switches( find(switches==1) );
     switches = (switches ==1);
    
    temp = find(switches);

    switches = switches';
    
    R = times(temp)';  

    R(:,2) = R(:,1)+winSize;
    
    result{i} = R;
end


%% Plots of output


figure(100);
plot(times,result{1})
