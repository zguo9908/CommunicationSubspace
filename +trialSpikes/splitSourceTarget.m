function [X_source, X_target, nSource, nTarget, index_source, index_target] ...
    = splitSourceTarget(sourceArea,numResult, X_hpc, X_pfc, varargin)
% Based on the number of directions asked, this function spilt the two
% firing matrices into source and target matices cell.

% INPUT
% numDirections: number of directions asked for (currently, can be set to 2
% or 4) e.g: hpc-hpc and hpc-pfc count as two
% numResult: number of resulting firing patterns
% X_region1/2: the firing matrices


% OUTPUT
% X_source: a cell of firing matrices that would later serve to predict the
% firing activities in the target regions.
% X_target: a cell of firing matrices that serve as the target

ip = inputParser;
ip.addParameter('specifiedSource', "CA1", @isstring);
ip.parse(varargin{:});
opt = ip.Results;

notUsingPFC = false;
notUsingHPC = false;

X_source = cell(1,numResult);
X_target = cell(2,numResult);
index_source = cell(1,numResult);
index_target = cell(2,numResult);

% number of cells in the two regions
[nhpc,~] = size(X_hpc{1});
[npfc,~] = size(X_pfc{1});



if npfc < 15 % not enough to be used as source
    notUsingPFC = true;
end

if nhpc < 15 % same above
    notUsingHPC = true;
end

if nhpc > npfc
    
    if notUsingHPC
        warning ("both regions have too few cells")
        nSourceHPC = -1;
        nSourcePFC = nSourceHPC;
    end
    
    if nhpc - npfc <15
        nSourceHPC = 15; % least # in source 
    else
        nSourceHPC = nhpc - npfc;
    end
    
    if notUsingPFC
        warning ("PFC has too few cells")
        nSourcePFC = -1;
    else
        nSourcePFC = ceil(0.5*npfc); 
    end
    
else
    
    if notUsingPFC
        warning ("both regions have too few cells")
        nSourceHPC = -1;
        nSourcePFC = nSourceHPC;
    end
    if npfc - nhpc <15
        nSourcePFC = 15;
    else
        nSourcePFC = npfc - nhpc;
    end
    
    if notUsingHPC
        warning ("HPC has too few cells")
        nSourceHPC = -1;

    else
        nSourceHPC = ceil(0.5*nhpc);
    end

end

if opt.specifiedSource == "CA1"

    nSource = nSourceHPC;
    if nSource == -1
        warning("not enough for source, maybe another animal");
    end
    nTarget = nhpc-nSource;
    for i = 1:numResult
        
        location_source = randperm(nhpc);
        location_source = location_source(1:nSource);
        location_source = sort(location_source);
        curr_source = X_hpc{i}(location_source,:);

        location_target1 = setdiff(1:nhpc, location_source);
     
        X_source{i} = curr_source;
        index_source{i} = location_source;
        index_target{1,i} = location_target1;
        
        curr_target =  X_hpc{i}(location_target1,:);
        X_target{1,i} = curr_target;
        
        
        location_target2 = randperm(npfc);
        location_target2 = location_target2(1:nTarget);
        curr_pfc = X_pfc{i}(location_target2,:);
        index_target{2,i} = location_target2;
        X_target{2,i} = curr_pfc(1:nTarget,:);
    end
else

    nSource = nSourcePFC;
    if nSource == -1
        warning("not enough for source, maybe another animal");
    end
    nTarget = npfc-nSource;
    for i = 1:numResult
          
        location_source = randperm(npfc);
        location_source = location_source(1:nSource);
        location_source = sort(location_source);
        curr_source = X_pfc{i}(location_source,:);

        location_target1 = setdiff(1:npfc, location_source);
     
        X_source{i} = curr_source;
        index_source{i} = location_source;
        index_target{1,i} = location_target1;
        
        curr_target =  X_pfc{i}(location_target1,:);
        X_target{1,i} = curr_target;
        
        
        location_target2 = randperm(nhpc);
        location_target2 = location_target2(1:nTarget);
        curr_hpc = X_hpc{i}(location_target2,:);
        index_target{2,i} = location_target2;
        X_target{2,i} = curr_hpc(1:nTarget,:);

    end
    
    
end
