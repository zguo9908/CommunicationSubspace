function [patterns, controls] = splitPatternsVSControls(all_structs, numPatterns, numControls,linearized)

% this function splits the structs array into one that's all actual
% patterns and another one that's all control patterns, assuming the last
% dimension of the sturct array represents the number of patterns and
% controls added

if ~linearized % when the Patterns struct array is first computed
    patterns = all_structs(:,:,numPatterns);
    controls = all_structs(:,:,numControls);
else
    numStructs = numel(all_structs);
    disp(numStructs)
    disp(numStructs * (numPatterns/(numPatterns+numControls)))
    patterns = all_structs(1:numStructs * (numPatterns/(numPatterns+numControls)));
    controls = all_structs(numStructs * (numPatterns/(numPatterns+numControls))+1: end);
end

