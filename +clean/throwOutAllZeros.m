function [firingmatrix, allZeroNeurons] = throwOutAllZeros(firingmatrix)

% this function takes in a firing rate matrix and throws out the 

allZeroNeurons = find(all(firingmatrix == 0,1));
firingmatrix(:,allZeroNeurons) = [];

end

