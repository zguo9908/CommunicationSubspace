function [timeAxis, times_spiking, spikeCountMatrix, spikeRateMatrix, areaPerNeuron] = getSpikeTrainSlow(animal, timebinSize, samplingRate)
%GETSPIKES Summary of this function goes here
%
% Input
% -----
% animal: name of the animal
%
% timebinSize : float
%   Size of time bins in seconds for spikeCountMatrix and
%   spikeRateMatrix. if given this option alone this determines both
%   sampling rate and window size.
%
% samplingRate : float
%   if given this (in addition to timebinSize), a time bin size window
%   is swept along the time axis collecting samples at the
%   samplingRate.
%
% Output
% -----
% times_spiking : cell
%   For every neuron, one list of spike times across the sessions,
%   stored in a cell
%
% spikeCountMatrix : timebins x neurons matrix
%   Spike counts per timebin
%
% spikeRateMatrix : timebins x neurons matrix
%   Spike rate per timebin
%
% areaPerNeuron : 1 x neurons string "" matrix
%   List of area per neuron. Used to later separate the columns
%   corresponding to each brain area. This function gathers that information
%   from the cellinfo file associated with the filename, by exacting animal
%   name from filename and appending 'cellinfo.mat'
%end

%% setting up inputs
%     filename = "JS15spikes01.mat";
%     timebinSize = 0.1;

load(animal + "spikes01.mat");
% % In this case using for-looping integers is easier than raw forlooping the
% % underlying structures. Easier than forlooping the underlying structures
session = 1;

%------------------------------------------------------------------
% to be implemented:
% given an animal get to the directory and get the spikes structure
%     animal = append(extractBetween(filename, 1, 4),"cellinfo.mat");
load (animal + "cellinfo.mat");

%------------------------------------------------------------------

spikes = spikes{session}; % single day epoches
cellinfo = cellinfo{session};
num_tetrode = numel(spikes{1});

cell_index = [];

for epoch = 1:numel(spikes)
    for tetrode = 1:numel(spikes{epoch})
        epoch_firing = spikes{epoch}{tetrode};
        for neuron = 1:numel(epoch_firing)
            if (numel(spikes{epoch}{tetrode}{neuron})~=0)
                cell_index(end+1,:) = [tetrode,neuron];
            end
        end
    end
end

% select unique cells
cell_index = unique(cell_index,"rows");

size_cellindex = size(cell_index);
num_cells = size_cellindex(1);


areaPerNeuron = string.empty;
times_spiking = cell(1,num_cells);


for epoch = 1:numel(spikes)
    for tetrode = 1:numel(spikes{epoch})
        epoch_firing = spikes{epoch}{tetrode};
        for neuron = 1:numel(epoch_firing)
            if (numel(spikes{epoch}{tetrode}{neuron})~=0)
                for i = 1:size_cellindex(1)
                    try
                        if neuron == cell_index(i,2) && tetrode == cell_index (i,1)
                            if isempty(spikes{epoch}{tetrode}{neuron}.data)
                                continue
                            else
                                neuron_data = spikes{epoch}{tetrode}{neuron}.data(:,1);
                                times_spiking{i} = [times_spiking{i}, neuron_data'];
                                
                                areaPerNeuron(i) = cellinfo{epoch}{tetrode}{neuron}.area;
                            end
                        end
                    catch
                        keyboard
                    end
                end
            end
        end
    end
end


%example figure of spiking time for one neuron
% figure (100)
% plot(times_spiking{1})
% xlabel("time of spiking")
%
% figure (101)
% plot(times_spiking{89})
% xlabel("time of spiking")

%% generate spike count/rate matrices
start_time = intmax;
end_time = -1;
for i = 1:num_cells
    try
        curr_start = times_spiking{i}(1);
    catch
        keyboard
    end
    curr_end = times_spiking{i}(end);
    
    if curr_start<start_time
        start_time = curr_start;
    end
    
    if curr_end>end_time
        end_time = curr_end;
    end
end

if nargin == 3 && ~isempty(samplingRate) %RY added new option for sampling rate + window size, copy-pasting some of your code and modifying
    
    samplingPeriod = 1/samplingRate;
    timeAxis = start_time:samplingPeriod:end_time;
    nTimes = length(timeAxis);
    spikeCountMatrix = zeros(num_cells,nTimes);
    
    
    chunkSize = 500;
    
    %         spikeCountSlice = zeros(1, nTime);
    %         for i = progress(1:num_cells,'Title','cells') % surrounding a 1:num_cells with progress() adds a progress bar
    %             for t = progress(1:chunkSize:numel(timeAxis),'Title','TimeChunks')
    %                 windowStarts = (timeAxis(t:min(t+chunkSize, nTime)) - timebinSize/2);
    %                 windowStops =  (timeAxis(t:min(t+chunkSize, nTime)) + timebinSize/2);
    %                 spikeCountSlice(t:min(t+chunkSize,nTime)) = ...
    %                     sum( times_spiking{i}' >= windowStarts & ...
    %                     times_spiking{i}' < windowStops);
    %             end
    %             spikeCountMatrix(i, :) = spikeCountSlice;
    %         end
    %         spikeRateMatrix = spikeCountMatrix/timebinSize;
    
    
    windowStarts = (timeAxis - timebinSize/2);
    windowStops =  (timeAxis + timebinSize/2);
    p = ProgressBar(num_cells, ...
        'Title', 'cells');
    %'IsParallel',true);
    p.setup([],[],[]);
    cleanupFunction = onCleanup(@(x) ProgressBar.deleteAllTimers());
    for i = 1:num_cells % surrounding a 1:num_cells with progress() adds a progress bar
        one_cell_spiking = times_spiking{i};
        nSpikes = numel(one_cell_spiking);
        spikeCountSlice = zeros(1, nTimes);
        for t = progress(1:chunkSize:numel(times_spiking{i}) ,'Title','SpikeChunks')
            spikeChunk = one_cell_spiking(t:min(t+chunkSize, nSpikes));
            spikeCountSlice = spikeCountSlice + ...
                sum( spikeChunk' >= windowStarts & ...
                spikeChunk' < windowStops);
        end
        spikeCountMatrix(i, :) = spikeCountSlice;
        p.step([], [], []);
    end
    spikeRateMatrix = spikeCountMatrix/timebinSize;
    
else % standard option, just window size alone
    
    timeAxis = start_time:timebinSize:end_time;
    
    
    spikeCountMatrix = zeros(num_cells,length(timeAxis)-1);
    spikeRateMatrix = zeros(num_cells,length(timeAxis)-1);
    
    for i = 1:num_cells
        [spike_count,timeAxis] = histcounts(times_spiking{i},timeAxis);
        spikeCountMatrix(i,:) = spike_count;
        spikeRateMatrix(i,:) = spike_count/timebinSize;
    end
end
end
