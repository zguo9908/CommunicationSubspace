function [X_regionPattern] = separateSpikes(spikeSampleMatrix, areaPerNeuron, region)
[~,nPatterns] = size(spikeSampleMatrix);
X_regionPattern = cell(1,nPatterns);
for r = 1:nPatterns
    curr = spikeSampleMatrix{r};
    [nNeurons,~] = size(curr);
    
    currResult = [];
    for i = 1:nNeurons

        if areaPerNeuron(i) == region
       
            currResult = [currResult;curr(i,:)];
        end
    end
    X_regionPattern{r} = currResult;
end


end