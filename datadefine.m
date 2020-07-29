function folder = datadefine()
% Automates pointing to the data  folders location on our computers

if ispc
    disp('Ziyi: insert your data folder path here');
    folder = '';
elseif isunix
    folder = '/Volumes/GoogleDrive/Share\ Drives/CommSubspace';
end