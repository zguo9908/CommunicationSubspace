function [X_source, X_target, nSource] = splitSourceTarget(numDirections, numResult, X_hpc, X_pfc, varargin)
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
ip.addParameter('specifiedSource', "HPC", @isstring);
ip.parse(varargin{:});
opt = ip.Results;

X_source = cell(numDirections/2,numResult);
X_target = cell(numDirections,numResult);

% number of cells in the two regions
[n1,~] = size(X_hpc{1});
[n2,~] = size(X_pfc{1});

if numDirections == 2 && opt.specifiedSource == "HPC"
    nSource = min(n1, abs(n1-n2));
end 
% elseif numDirections ==2 && opt.specifiedSource == "PFC"
%     nSource = min(n1, abs(n1-n2));
% else
%     nSource1 = 
for i = 1:numResult
    curr_hpc = X_hpc{i}(randperm(n1),:);
    curr_source = curr_hpc(1:nSource,:);
    X_source{i} = curr_source;
    curr_target = curr_hpc(nSource+1:max(n1, n2),:);
    X_target{1,i} = curr_target;
    X_target{2,i} = X_pfc{i};
end

end

