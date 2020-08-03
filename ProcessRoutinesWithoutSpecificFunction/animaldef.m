function animalinfo = animaldef(animalname)

% Get the user name! and determine user specific root folder
if ispc
    rootfolder = "\\citadel.bio.brandeis.edu\sharespace-commsub\";
elseif  ismac
    rootfolder = "~/Data/commsubspace/";
end

switch animalname
    
    % Goal Maze Animals
    case 'RY7', animalinfo={'RY7',['/media/ryoung/Thalamus' filesep 'ry_GoalCoding_Project/RY7_experiment/RY7_direct'],'RY7'};
    case 'RY9', animalinfo={'RY9',['/media/ryoung/Thalamus' filesep 'ry_GoalCoding_Project/RY9_experiment/RY9_direct'],'RY9'};
        
        % Phase disrupt animals
    case 'CL1', animalinfo={'CL1',[rootfolder 'Phase/CL1_direct/'],'CL1'};
    case 'SG7', animalinfo={'SG7',[rootfolder 'OdorPlace/SG7Expt/SG7_direct/'],'SG7'};
        
        % CA3-CA1 animals
        % Right now, just animals we navigated tetrodes down to CA3 to figure out how to get there
    case 'EG1'
        animalinfo = {'EG1', [rootfolder 'CA3-CA1/Efizz/EG1_direct/'] , 'EG1'};
        
        
        % Single day animals
    case 'YD6'
        animalinfo = {'YD6', strcat(rootfolder, 'W-Track_WellTrained_EPHYS/YD6_direct/'),'YD6'};
    case 'ER1'
        animalinfo = {'ER1', strcat(rootfolder, 'SingleDayExpt/ER1_NEW_direct/'), 'ER1'};
    case 'KL8'
        animalinfo = {'KL8', strcat(rootfolder, 'SingleDayExpt/KL8_direct/'), 'KL8'};
    case 'JS7'
        animalinfo = {'JS7', strcat(rootfolder, 'W-track_SingleDay/Efizz/JS7_direct/'),'JS7'};
    case 'SJ6'
        animalinfo = {'SJ6', strcat(rootfolder, 'W-track_SingleDay/Efizz/SJ6_direct/'),'SJ6'};
    case 'SJ1'
        animalinfo = {'SJ1', strcat(rootfolder, 'W-track_SingleDay/Beh/SJ1_direct/'),'SJ1'};
    case 'SJ2'
        animalinfo = {'SJ2', strcat(rootfolder, 'W-track_SingleDay/Beh/SJ2_direct/'),'SJ2'};
    case 'SJ3'
        animalinfo = {'SJ3', strcat(rootfolder, 'W-track_SingleDay/Beh/SJ3_direct/'),'SJ3'};
    case 'SJ4'
        animalinfo = {'SJ4', strcat(rootfolder,  'W-track_SingleDay/Beh/SJ4_direct/'),'SJ4'};
    case 'JS12'
        animalinfo = {'JS12', strcat(rootfolder, 'SingleDayExpt/JS12_direct/'),'JS12'};
    case 'JS13'
        animalinfo = {'JS13', strcat(rootfolder, 'SingleDayExpt/JS13_direct/'),'JS13'};
    case 'JS14'
        animalinfo = {'JS14', strcat(rootfolder, 'SingleDayExpt/JS14_direct/'),'JS14'};
    case 'JS15'
        animalinfo = {'JS15', strcat(rootfolder, 'SingleDayExpt/JS15_direct/'),'JS15'};
    case 'JS17'
        animalinfo = {'JS17', strcat(rootfolder, 'SingleDayExpt/JS17_direct/'),'JS17'};
    case 'JS21'
        animalinfo = {'JS21', strcat(rootfolder, 'SingleDayExpt/JS21_direct/'),'JS21'};
    case 'JS21'
        animalinfo = {'JS21', strcat(rootfolder, 'SingleDayExpt/JS21_direct/'),'JS21'};
    case 'ZT2'
        animalinfo = {'ZT2', strcat(rootfolder, 'SingleDayExpt/ZT2_direct/'),'ZT2'};
        
        
        % New 8day HC-PFC animals
    case 'SJ5'
        animalinfo = {'SJ5', strcat(rootfolder, 'HP_8dayExpt/SJ5_direct/'),'SJ5'};
        
        % Ripple interruption, old animals
    case 'sjc'
        animalinfo = {'sjc', strcat(rootfolder, 'RippleDisruption_all/RippleDisruption/sjc_direct/'), 'sjc'};
    case 'RE1'
        animalinfo = {'RE1', strcat(rootfolder, 'RippleDisruption_all/RippleDisruption/RE1_direct/'), 'RE1'};
    case 'RNa'
        animalinfo = {'RNa', strcat(rootfolder, 'RippleDisruption_all/RippleDisruption/RNa_direct/'), 'RNa'};
    case 'RNb'
        animalinfo = {'RNb', strcat(rootfolder, 'RippleDisruption_all/RippleDisruption/RNb_direct/'), 'RNb'};
    case 'RNc'
        animalinfo = {'RNc', strcat(rootfolder, 'RippleDisruption_all/RippleDisruption/RNc_direct/'), 'RNc'};
    case 'RNd'
        animalinfo = {'RNd', strcat(rootfolder, 'RippleDisruption_all/RippleDisruption/RNd_direct/'), 'RNd'};
    case 'RCa'
        animalinfo = {'RCa', strcat(rootfolder, 'RippleDisruption_all/RippleDisruption/RCa_direct/'), 'RCa'};
    case 'RCb'
        animalinfo = {'RCb', strcat(rootfolder, 'RippleDisruption_all/RippleDisruption/RCb_direct/'), 'RCb'};
    case 'RCc'
        animalinfo = {'RCc', strcat(rootfolder, 'RippleDisruption_all/RippleDisruption/RCc_direct/'), 'RCc'};
    case 'RCd'
        animalinfo = {'RCd', strcat(rootfolder, 'RippleDisruption_all/RippleDisruption/RCd_direct/'), 'RCd'};
    case 'REc'
        animalinfo = {'REc', strcat(rootfolder, 'RippleDisruption_all/RippleDisruption/REc_direct/'), 'REc'};
    case 'REd'
        animalinfo = {'REd', strcat(rootfolder, 'RippleDisruption_all/RippleDisruption/REd_direct/'), 'REd'};
    case 'REe'
        animalinfo = {'REe', strcat(rootfolder, 'RippleDisruption_all/RippleDisruption/REe_direct/'), 'REe'};
    case 'REf'
        animalinfo = {'REf', strcat(rootfolder, 'RippleDisruption_all/RippleDisruption/REf_direct/'), 'REf'};
    case 'REg'
        animalinfo = {'REg', strcat(rootfolder, 'RippleDisruption_all/RippleDisruption/REg_direct/'), 'REg'};
    case 'REh'
        animalinfo = {'REh', strcat(rootfolder, 'RippleDisruption_all/RippleDisruption/REh_direct/'), 'REh'};
        
        % Hippocampal-prefrontal animals
    case 'HPa'
        animalinfo = {'HPa', strcat(rootfolder, 'HP_8dayExpt/HPExpt/HPa_direct/'), 'HPa'};
    case 'HPb'
        animalinfo = {'HPb', strcat(rootfolder, 'HP_8dayExpt/HPExpt/HPb_direct/'), 'HPb'};
    case 'HPc'
        animalinfo = {'HPc', strcat(rootfolder, 'HP_8dayExpt/HPExpt/HPc_direct/'), 'HPc'};
    case 'Nadal'
        animalinfo = {'Ndl', strcat(rootfolder, 'HP_8dayExpt/HPExpt/Ndl_direct/'), 'Ndl'};
    case 'Rosenthal'
        animalinfo = {'Rtl', strcat(rootfolder, 'HP_8dayExpt/HPExpt/Rtl_direct/'), 'Rtl'};
    case 'Borg'
        animalinfo = {'Brg', strcat(rootfolder, 'HP_8dayExpt/HPExpt/Brg_direct/'), 'Brg'};
        
        % UH  OH -- animal name not recognized ...
    otherwise
        error(['Animal ',animalname, ' not defined.']);
end

animalinfo{2} = char(animalinfo{2});

if ~exist(animalinfo{2},'dir')
    warning(['Directory ' animalinfo{2} ' not found! Creating...']);
    mkdir(animalinfo{2});
end

fprintf('Using animal info at %s \n', animalinfo{2});
