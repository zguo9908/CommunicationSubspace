function folder = codedefine()

if ispc
    disp('Defining code folder for Ziyi machine');
    folder = 'C:\\Users\BrainMaker\Matlab Drive\Shared\';
elseif ismac
    disp('Defining code folder for Ryan''s machine');
    folder = '~/Data/Matlab-DRIVE/Shared/';
elseif isunix
    disp('Defining code folder for Ryan''s linux machine');
    folder = '/Volumes/';
end
