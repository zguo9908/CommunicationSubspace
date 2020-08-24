clear
%% Paths
addpath(genpath(pwd)) % all folders in utils added, including semedo code
if ispc
    paths = " C:\Users\BrainMaker\commsubspace\SingleDayExpt";
elseif ismac
    paths(1) = "/Volumes/sharespace-commsub/data";
    paths(2) = "~/Data/commsubspace";
end

arrayfun(@(path) addpath(genpath(path)), paths);
% %%
% FIELD = [];
% numToRun = numel(FIELD);
%% Script parameters
% -----------------------------------------------------------------
Option = struct();

% Option.animal = [];
% Option.generateFromRipTimes = true; 
% Option.generateH = "fromFilteredEEG "+" fromRipTimes";
% Option.samplingRate  = "empty" ;         % For spikes.getSpikeTrain
% Option.spikeBinSize  = 0.1;          % 100 milliseconds
% Option.timesPerTrial = 10;         % 10 times per trial
% Option.winSize       = {[-0.15, 0.15]};              % size of the window
% Option.sourceArea    = "PFC";
% Option.equalWindowsAcrossPatterns = true;    % whether all three patterns have the same #windows
% Option.singleControl = true;                 % whether to use just one control column
% Option.numPartition = 10;                    % ways to split source and target
% usingSingleprediction = true;
% 
% winSize = Option.winSize{1};

animals = ["JS15"];
% animals = ["JS13", "JS14","JS21", "ER1", "KL8"];

for iAnimal = 1:numel(animals)
    Option.animal = animals(iAnimal)
    TheScript
end 