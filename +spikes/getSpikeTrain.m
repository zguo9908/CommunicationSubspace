function [timeAxis, times_spiking, spikeCountMatrix, spikeRateMatrix] = getSpikeTrain(filename,timebinSize)
%GETSPIKES Summary of this function goes here
%
% Input
% -----
% animal : string
%   Name of the animal
%
% session : number
%   Number of the session (day)
%
% timebinSize : number
%   Size of time bins in seconds for spikeCountMatrix and
%   spikeRateMatrix.
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
%end

%% setting up inputs

load(filename);
% % In this case using for-looping integers is easier than raw forlooping the
% % underlying structures. Easier than forlooping the underlying structures
session = 1;

%------------------------------------------------------------------
% to be implemented:
% given an animal get to the directory and get the spikes structure
%------------------------------------------------------------------

spikes = spikes{session}; % single day epoches
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
times_spiking = cell(1,num_cells);
for epoch = 1:numel(spikes)
    for tetrode = 1:numel(spikes{epoch})
        epoch_firing = spikes{epoch}{tetrode};
        for neuron = 1:numel(epoch_firing)
           if (numel(spikes{epoch}{tetrode}{neuron})~=0)
                for i = 1:size_cellindex(1)
                    if neuron == cell_index(i,2) && tetrode == cell_index (i,1)
                        neuron_data = spikes{epoch}{tetrode}{neuron}.data(:,1);
                        times_spiking{i} = [times_spiking{i}, neuron_data'];
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
    curr_start = times_spiking{i}(1);
    curr_end = times_spiking{i}(end);
    
    if curr_start<start_time
        start_time = curr_start;
    end
    
    if curr_end>end_time
        end_time = curr_end;
    end
end
timeAxis = start_time:timebinSize:end_time;

spikeCountMatrix = zeros(num_cells,length(timeAxis)-1);
spikeRateMatrix = zeros(num_cells,length(timeAxis)-1);

for i = 1:num_cells
    [spike_count,timeAxis] = histcounts(times_spiking{i},timeAxis);
    spikeCountMatrix(i,:) = spike_count;
    spikeRateMatrix(i,:) = spike_count/timebinSize;
end
end