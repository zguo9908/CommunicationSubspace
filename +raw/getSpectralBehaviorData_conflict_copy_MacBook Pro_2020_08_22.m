function data = getSpectralBehaviorData(animalList, varargin)
% GETSPECTRALDATA not a particularly pretty function but fetches spectral and positional data
% in one large data struct. I wrote it for SeqNMF stuff, but have since adapted it into
% a function for more general use.
%
% Inputs
% ------
%
% animalList : cell of char or array of str
%   If you ask for more than one animal, it will attempt to
%   change the times so that they're unique per animal.
% 
% (Optional key-word args)
% -------------------------
%
% epochType  : str, default='run', possibilities={all|run|sleep}
%
% debug : logical, default=false
%
% subsample : numeric, default=false
%   Whether to subsample random chunks of each animals epoch. If
%   subsample is a number 0 to 1, its the fraction of an epoch to
%   randomly, but contiguosly sample.
%
% storeDatBeforeSubsamp : logical, default=false
%   If subsampling, should we store away the original data?
%
% fields : cell of str
%   Which data fields to sample per aninmal
%
%       Efizz fields:
%       'C'
%       'S1'
%       'S2'
%       'wpli'
%       'phi'
%       'phi_sin'
%       'phi_cos'
%
%       Behavior fields:
%       'velocity'
%       'trajdist'
%       'trajbound'
%       'reward'
%
%
% Output
% ------
%
% data

% --------------------------------
% Parse optional keyword-arg pairs
% --------------------------------

isNew = false;
ip = inputParser;
ip.addParameter('epochType', 'run', @(x) ischar(x) || isstring(x))
ip.addParameter('fields', {'C','wpli','S1','S2'}); 
ip.addParameter('debug', false)
ip.addParameter('subsample', false); 
ip.addParameter('removeInterepochPeriod', false);
ip.addParameter('removeInteranimalPeriod', false);
ip.addParameter('storeDatBeforeSubsamp', false);
ip.addParameter('onehotencoding', false);
ip.addParameter('fieldpack_kws', {});
ip.parse(varargin{:});
opt = ip.Results;

%% Loop over and collect efizz/behavior data

P = []; % for aggreagting pos
C = []; % for aggreagting efizz
T = []; % for aggregating lregion

cnt = 0;
animal_cnt = 0;
offset = 0;
dayepoch_prev = [0 0];

for animal = animalList

    animal_cnt = animal_cnt + 1;

    % Determine the days and epochs for the animal
    animaldef_ = animaldef(animal{1});
    animalToPath(animal{1});

    task = loaddatastruct(animaldef_{2:3}, 'task');
    search_results = cellfetch(task, 'environment');

    switch opt.epochType
    case 'run'
        epoch_filter = cellfun(@(x) isequal(x,'wtr1'), search_results.values);
    case 'sleep'
        epoch_filter = cellfun(@(x) isempty(x), search_results.values);
    case 'all'
        epoch_filter = ones(1,numel(search_results.values));
    otherwise
        error("Improper epoch type")
    end
    dayepochs = search_results.index(epoch_filter,:);
    if opt.debug
        dayepochs
    end

    % ------------------------------
    % Load up all data for an animal
    % ------------------------------
    for dayepoch = dayepochs'

        day   = dayepoch(1);
        epoch = dayepoch(2);
        if sum(dayepoch - dayepoch_prev) < 0
            load(prev_avgeeg_file)
            avgeeg = ffend(avgeeg);
            if opt.debug
                offset = offset  + avgeeg.endtime
            end
        end
        % ------------------------
        % -- load spectral info --
        % ------------------------
        if isequal(opt.epochType,'run')
            if isNew
                cgramfile = "cgramcnew"
            else
                cgramfile = 'cgramc';
            end
        else
           if isNew
                cgramfile = "cgramcnew"
           else
                cgramfile = 'cgramc';
            end

        end
        cgramfile = "cgramc";
        temp = join([animal{1} cgramfile sprintf('-%02d-%02d.mat', dayepoch)]);
        cgramc_file = strrep(temp, ' ', '');
        if ~exist(cgramc_file, 'file')
            disp('File does not exist .. Continuing')
            continue
        else
            disp("Processing")
        end
        load(cgramc_file)
        if isNew
        cgramc = cgramcnew;
        end
        data = ffend(cgramc);
        data.orig.t    = data.t;
        data.t         = data.t + offset;
        data.animalcnt = animal_cnt * ones(size(data.t));
        data.epoch     = epoch * ones(size(data.t));
        %data.t = data.t - avgeeg.starttime;
        
        % ------------------------
        % -- load position info --
        % ------------------------
        try
            C = [C; data];
            if ~isequal(opt.epochType, 'sleep')
                load([animal{1} 'pos' sprintf('%02d', day)])
                load([animal{1} 'lregion' sprintf('%02d', day)])
                p = pos{day}{epoch};
                p.data(:,1) = p.data(:,1) + offset;
                P = [P; repmat([animal_cnt, cnt],length(p.data),1), p.data];
                p = lregion{day}{epoch};
                p.time_orig = p.time; 
                p.time = p.time + offset;
                % Order of variables
                % animal epoch time trial region rewarded trajbound lindist trajdist
                T = [T; ...
                    repmat([animal_cnt,     cnt], length(p.time),1), p.time, ...
                    p.TR, p.TB, p.lindist, p.lindist .* sign(p.TB(:,2)-0.5), ...
                    p.tperf_timewise, p.time_orig]; 
            end
        catch ME
            warning(sprintf('Skipping animal %s, day %d epoch %d\n', animal{1}, dayepoch));
        end
        %L =  round(timescale/data.params.movingwin(2));
        prev_avgeeg_file = [animal{1} 'avgeeg' sprintf('%02d-%02d.mat', dayepoch)];
        dayepoch_prev = dayepoch;
    end
end

%% Store general
% Store potentially modified field list
clear data;
data.fields = opt.fields;


%% Store behavior
if isequal(opt.epochType,'run') && ~(isempty(P) || isempty(T)) 
    % -----------------------------------------
    % Convert pos and trial to common time
    % -----------------------------------------
    throwaway = ~ismember(P(:, 1:3), T(:, 1:3), 'rows');
    P(throwaway,:) = [];
    throwaway = ~ismember(T(:, 1:3), P(:, 1:3), 'rows');
    T(throwaway,:) = [];
    % Sort by animal epoch time
    P = sortrows(P, [1,2,3]);
    T = sortrows(T, [1,2,3]);
    % -------------------------------
    % Create convenience table arrays
    % -------------------------------
    assert(isequal(P(:,1:3),T(:,1:3)))
    % Behavior dataset
    P = num2cell(P,1);
    P = table(P{:}, ...
    'VariableNames', {'animal','epoch','time','x','y','dir','vel'});
    T = num2cell(T(:,4:end),1);
    T = table(T{:}, ...
    'VariableNames', {'trial','region','rew','traj', 'lindist', 'trajdist', 'tperf', 'orig_time'});
    behavior = [P,T];
    %data.T = [0 diff(data.t)];
    %data.T(data.T<0) = 0;
    %data.T = cumsum(data.T);
    %behavior.time(:) = [0; diff(behavior.time(:))];
    %behavior.time(behavior.time<0) = 0;
    %behavior.time = cumsum(behavior.time);
    data.behavior = behavior;
end



%% Electrophysiology fields
% Time and frequency
assert(~isempty(C),'C is empty')
data.efizz.t     = cat(2, C.t);
torig      = arrayfun(@(x) x.orig.t, C, 'UniformOutput', false);
data.efizz.f     = C.f;
clear torig;
% Animal
data.efizz.animalcnt = cat(2, C.animalcnt);
data.efizz.animalnames = string(animalList);
data.efizz.epoch     = cat(2, C.epoch);
% Other fields
for field = opt.fields
    if isfield(C,field{1})
        data.efizz.(field{1}) = cat(1, C.(field{1}));
    end
end
if any(contains(opt.fields, 'phi'))
    data.efizz.phi_sin = sin(data.phi);
    data.efizz.phi_cos = cos(data.phi);
    data.fields{cellfun(@(x) isequal(x, 'phi'), data.fields)} = 'phi_sin';
    data.fields = [data.fields 'phi_cos'];
end


%% Other data post-processing options
% ----------------------------------------
% Compute total seconds of data
% ----------------------------------------
Dt        = diff(data.efizz.t);
if opt.removeInteranimalPeriod
    Dt(Dt<0)  = 0;
end
if opt.removeInterepochPeriod
    Dt(Dt>10) = 0;
end
disp(['Total seconds of data: ' num2str(sum(Dt))])
clear Dt

% ------------------------------
% ADD HOTENCODED BEHAVIOR FIELDS
% ------------------------------
if opt.onehotencoding
    if any(contains(opt.fields, 'velocity'))
        behavior2data_inds = lookup(double(data.t), double(behavior.time));
        data.hoteencoding.behavior2data_inds = behavior2data_inds;
        data.hoteencoding.velocity = log(behavior.vel(behavior2data_inds));
        data.hoteencoding.velocity_repelem = 3;
        data.hoteencoding.orig_velocity    = data.hotencoding.velocity;
        data.hoteencoding.velocity         = label2mat(data.hotencoding.velocity, 20, 'method', 'quantile')';
        data.hoteencoding.velocity         = repelem(data.hotencoding.velocity, 1, data.hotencoding.velocity_repelem);
    end
    if any(contains(opt.fields, 'trajdist'))
        behavior2data_inds = lookup(double(data.t), double(behavior.time));
        data.hoteencoding.behavior2data_inds = behavior2data_inds;
        data.hoteencoding.trajdist = behavior.trajdist(behavior2data_inds);
        data.hoteencoding.trajdist_repelem = 3;
        data.hoteencoding.orig_trajdist    = data.hotencoding.trajdist;
        data.hoteencoding.trajdist         = label2mat(data.hotencoding.trajdist, 40)';
        data.hoteencoding.trajdist         = repelem(data.hotencoding.trajdist, 1, data.hotencoding.trajdist_repelem);
    end
    if any(contains(opt.fields, 'lindist'))
        behavior2data_inds = lookup(double(data.t), double(behavior.time));
        data.hotencoding.lindist = behavior.lindist(behavior2data_inds);
        data.hotencoding.lindist_repelem = 3;
        data.hotencoding.orig_lindist    = data.hotencoding.lindist;
        data.hotencoding.lindist         = label2mat(data.hotencoding.lindist, 20)';
        data.hotencoding.lindist         = repelem(data.hotencoding.lindist, 1, data.hotencoding.lindist_repelem);
    end
    if any(contains(opt.fields, 'aglindist'))
        behavior2data_inds = lookup(double(data.t), double(behavior.time));
        data.hotencoding.behavior2data_inds = behavior2data_inds;
        data.hotencoding.aglindist = behavior.lindist(behavior2data_inds);
        inbound = behavior.traj(behavior2data_inds) == 1;
        data.hotencoding.aglindist(inbound) = -1*data.hotencoding.aglindist(inbound) + 1; % Map 1->0 traj to 0->1, like the outbound
        data.hotencoding.aglindist_repelem = 3;
        data.hotencoding.orig_aglindist    = data.hotencoding.aglindist;
        data.hotencoding.aglindist         = label2mat(data.hotencoding.aglindist, 10)';
        data.hotencoding.aglindist         = repelem(data.hotencoding.aglindist, 1, data.hotencoding.aglindist_repelem);
       
    end
    if any(contains(opt.fields, 'trajbound'))
        behavior2data_inds = lookup(double(data.t), double(behavior.time));
        data.hotencoding.behavior2data_inds = behavior2data_inds;
        [~, ~, behavior] = characterize.futurePastTrajbound(behavior, 'n', 2);
        data.hotencoding.trajbound = [behavior.pasttrajbound(behavior2data_inds, :), behavior.rew(behavior2data_inds), behavior.futuretrajbound(behavior2data_inds, :)];
        if onehotencoding
            data.hotencoding.trajbound_repelem = 1;
            data.hotencoding.orig_trajbound    = data.hotencoding.trajbound;
            %data.trajbound         = repelem(data.trajbound, 1, data.trajbound_repelem);
            data.hotencoding.trajbound         = arrayfun(@(col) label2mat(data.hotencoding.trajbound(:,col)+1',2,'method','raw'), ...
                                    1:size(data.hotencoding.trajbound,2), 'UniformOutput', false);
            data.trajbound         = cat(1, data.trajbound{:})';
        end
    end
    if any(contains(opt.fields, 'reward'))
        behavior2data_inds = lookup(double(data.t), double(behavior.time));
        data.hotencoding.behavior2data_inds = behavior2data_inds;
        [~, ~, behavior] = characterize.futurePastReward(behavior, 'n', 2);
        data.hotencoding.reward = [behavior.pastreward(behavior2data_inds,:), behavior.rew(behavior2data_inds), behavior.futurereward(behavior2data_inds,:)];
        data.hotencoding.reward_repelem = 1;
        data.hotencoding.orig_reward    = data.hotencoding.reward;
        %data.reward         = repelem(data.reward, 1, data.reward_repelem);
        data.hotencoding.reward         = arrayfun(@(col) label2mat(data.hotencoding.reward(:,col)+1',2,'method','raw'), ...
                                1:size(data.hotencoding.reward,2), 'UniformOutput', false);
        data.hotencoding.reward         = cat(1, data.hotencoding.reward{:})';
    end

end
% -----------
% Pack fields
% -----------
if ~isempty(opt.fieldpack_kws)
    fieldpack_tmp = [opt.fieldpack_kws, 'groups', data.animalcnt];
    data.data = seq.packfields(data, data.fields, fieldpack_tmp{:})';
end

% -------------------------------------
% Store copy of data before subsampling
% -------------------------------------
if opt.storeDatBeforeSubsamp && opt.subsample > 0
    for field = union(data.fields, ...
                     {'S1','S2','plv','wpli','C','trajdist','f','t','data'})
        if ~ismember(field{1}, behavior.Properties.VariableNames)
            data.pre.(field{1}) = data.(field{1});
        end
    end
end

% ---------
% SUBSAMPLE
% ---------
if opt.subsample
    disp('Before')
    size(data.data)
    data = seq.subsampleByGroup(data, fields, opt.subsample, 'storeOrig', storeDatBeforeSubsamp);
    disp('After')
    size(data.data)
    assert(numel(data) > 0);
end


