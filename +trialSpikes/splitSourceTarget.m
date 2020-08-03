function [X_source, X_target, nSource] = splitSourceTarget(sourceArea,...
    numResult, X_hpc, X_pfc, varargin)
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

% number of cells in the two regions
[nhpc,~] = size(X_hpc{1});
[npfc,~] = size(X_pfc{1});

if npfc < 15
    notUsingPFC = true;
end

if nhpc < 15
    notUsingHPC = true;
end

if nhpc > npfc
    
    if notUsingHPC
        warning ("both regions have too few cells")
    end
    
    if nhpc - npfc <15
        nSourceHPC = 15;
    else
        nSourceHPC = nhpc - npfc;
    end
    
    if notUsingPFC
        warning ("PFC has too few cells")
    else
        nSourcePFC = floor(0.5*npfc); 
    end
else
    
    if notUsingPFC
        warning ("both regions have too few cells")
    end
    if npfc - nhpc <15
        nSourcePFC = 15;
    else
        nSourcePFC = npfc - nhpc;
    end
    
    if notUsingHPC
        warning ("HPC has too few cells")
    else
        nSourceHPC = floor(0.5*nhpc);
    end
    
end

if opt.specifiedSource == "CA1"
    nSource = nSourceHPC;
    for i = 1:numResult
        curr_hpc = X_hpc{i}(randperm(nhpc),:);
        curr_source = curr_hpc(1:nSource,:);
        X_source{i} = curr_source;
        curr_target = curr_hpc(nSource+1:nhpc,:);
        X_target{1,i} = curr_target;
        X_target{2,i} = X_pfc{i};
    end
else
    nSource = nSourcePFC;
    for i = 1:numResult
        curr_pfc = X_pfc{i}(randperm(npfc),:);
        curr_source = curr_pfc(1:nSource,:);
        X_source{i} = curr_source;
        curr_target = curr_pfc(nSource+1:npfc,:);
        X_target{1,i} = curr_target;
        X_target{2,i} = X_hpc{i};
    end
    
    
end
