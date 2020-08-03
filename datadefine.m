function folder = datadefine()
% Automates pointing to the data  folders location on our computers

if ispc
    disp('Defining data folder for Ziyi machine');
    folder = 'C:\\Users\BrainMaker\commsubspace\';
elseif ismac
    disp('Defining data folder for Ryan''s machine');
    folder = '~/Data/commsubspace/';
elseif isunix
    disp('Defining data folder for Ryan''s linux machine');
    folder = '/Volumes/';
end
