%% Parse optional arguments in varargin
% (See docfile on inputParser
%  https://www.mathworks.com/help/matlab/ref/inputparser.html#d120e629568)
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

switch parsedOptionalInput.controlType
    case 'none'
        % Do not add a control column
    case 'shuffle'
        H = generateControlColumn(H, times, spectrogram);
    case 'shuffle_avoidpatterns'
        H = generateControlColumn(H, possible_times, spectrogram);
    otherwise
        error("controlType has improper value = " + parsedOptionalInput.controlType)
end

%% shuffle
X(randperm(length(X)))

%% interpret more points based on few known

%using linespace to generate equidistance points as query points
