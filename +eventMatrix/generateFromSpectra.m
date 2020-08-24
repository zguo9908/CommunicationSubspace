function [H, Hvals, Hnanlocs, times] = generateFromSpectra(times, spectrogram, frequencyAxis, ...
                                 frequenciesPerPattern, varargin)

    % generates H matrix from lfp data
    
    % Inputs
    % ------
    % times: 1 x T matrix of all the possible time points
    %
    % spectrogram : T x F matrix
    %   Matrix recording the power of each frequency band at each time.
    %   Where T is the number of times and F number of frequencies.
    %
    % frequencyAxis : 1 x F
    %   Gives the actual frequency-axis, the values in Hz for each of the
    %   1:F frequncies.
    %
    % frequenciesPerPattern : [freqRange1; ...; freqRangeN] 2xN matrix
    %   2xN matrix of frequncy ranges for N different bands/patterns. These
    %   specify the frequency ranges for the columns of H.
    %
    % Output
    % ------
    % H : T x N matrix
    %   Average amount of power for the N patterns/bands for T times.
    % Hvals: same as H to keep things consistent with other generateH
    % functions
    % Hnanlocs: encoding where the nans are with nan and the others with 1
    
    ip = inputParser;
    ip.addParameter('controlType','none', ... % name of an optional input and it's default value
        @(x) ischar(x) || isstring(x)); % test run on optional input (if it fails this test, it triggers an error)
    ip.parse(varargin{:})
    parsedOptionalInput = ip.Results;
    
    % This code adds 'controlType' as an optional input. To call a funciton
    % with optional input, write:
    %
    % makeWindows(times, Q, winSize, 'controlType', 'shuffle')
    %
    % To run without an optional (it defaults to 'none'), write:
    %
    % makeWindows(times, Q, winSize)
    % 
    % This is same as calling:
    %
    % makeWindows(times, Q, winSize, 'controlType', 'none')
    
    %% setting up inputs
    
    % geting sizes
    [nTimes,~] = size(spectrogram);
    [nPattern,~] = size(frequenciesPerPattern);
    
    % instantiate H and possible_times
    H = zeros(nTimes,nPattern);
    possible_times = cell(1,nPattern);
    
    for i = 1:nPattern
        
        f1 = frequenciesPerPattern(i,1);
        f2 = frequenciesPerPattern(i,2);
        frequency_subset = frequencyAxis > f1 & frequencyAxis < f2;
        result = spectrogram(:,frequency_subset); %frequencies that fit the range
        % result is T by (number of discrete frequencies that fit in the
        % range)
        
        result = result'; % the average power of each discrete frequency, since mean operates on a row-wise fashion
        % we need to transpose 
        H(:,i) = (mean(result))';
        
        if any(isnan(H(:,i)))
            keyboard
        end
        curr_times = times(frequency_subset);
        possible_times{i} = setdiff(times,curr_times); % exclude times for the next pattern
    end      

Hvals = H;
Hnanlocs = ones(size(H));