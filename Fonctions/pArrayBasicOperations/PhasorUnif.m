function [paddedPhasors,sizesMatrix] = PhasorUnif(phasorArray)
% PHASORUNIF takes an arbitrary number of phasor arrays and pads their third
% dimension (representing harmonics) to match the maximum number of harmonics
% present in the input arrays. The padding is done with zeros and is added
% symmetrically to the beginning and end of the third dimension. If any of the
% input arrays are of type 'PhasorArray', their value property is used for
% padding and the result is wrapped back into a 'PhasorArray'. The function
% preserves 'sdpvar' and 'sym' types in the input arrays.
arguments (Repeating)
    phasorArray %accept an arbitrary number of phasor arrays as input
end

    numArrays=numel(phasorArray);

% Check if elements of phasorArray are of type 'PhasorArray' and replace them with their values
isPhasorArray=cellfun(@(x) isa(x,'PhasorArray'),phasorArray);
phasorArray(isPhasorArray)=cellfun(@(x) x.value,phasorArray(isPhasorArray),'UniformOutput',false);

% Get sizes of all arrays
sizesMatrix = zeros(numArrays, 3);
for i = 1:numArrays
    sizesMatrix(i, :) = [size(phasorArray{i}, 1), size(phasorArray{i}, 2), size(phasorArray{i}, 3)];
end

% Calculate maximum size in the third dimension
maxThirdDimSize=max(sizesMatrix(:,3));

% Calculate padding needed for each array to match max size in third dimension
targetSize=ones(size(sizesMatrix,1),1)*[sizesMatrix(1,1) sizesMatrix(1,2) maxThirdDimSize];
paddingSizes=(targetSize-sizesMatrix)/2;
paddingSizes=num2cell(paddingSizes,2)';

% Initialize output cell array
paddedPhasors=cell(numArrays,1)';

% Loop over each array in phasorArray
for i = 1:numel(phasorArray)
    % Get current array
    currentArray=phasorArray{i};
    
    % Create padding of zeros
    padding = zeros(size(currentArray,1), size(currentArray,2), paddingSizes{i}(3));
    
    % Concatenate padding, array, and padding along the third dimension
    paddedPhasors{i} = cat(3, padding, currentArray, padding);
end

% Replace original PhasorArray elements with padded ones
paddedPhasors(isPhasorArray)=cellfun(@(x) PhasorArray(x,reduce=false),paddedPhasors(isPhasorArray),'UniformOutput',false);
end