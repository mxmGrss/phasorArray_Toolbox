function [summedPhasors, paddedPhasors] = PhasorArrayAdd(inputPhasors)
    % PHASORARRAYADD performs padding (and centering) of an array of phasors and then
    % adds their slices to form A(t) + B(t) in the harmonic domain.
    % inputPhasors is a cell array of 3D arrays to be added.
    % The function returns a 3D array summedPhasors that is the sum of the input arrays,
    % and a cell array paddedPhasors that contains the input arrays after padding.
    arguments (Repeating)
        inputPhasors
    end
  
    % Uniformize the size of the input arrays and get their original sizes
    [paddedPhasors, originalSizes] = PhasorUnif(inputPhasors{:});
    
    % Check that all arrays have the same size in the first and second dimensions

    if range(originalSizes(:,1))~=0 || range(originalSizes(:,2))~=0
        error("Non compatible dimensions. All arrays must have the same size in the first and second dimensions.")
    end

    % Sum all arrays along the fourth dimension to form the output array
    summedPhasors=sum(cat(4,paddedPhasors{:}),4);
end

function r = range(X)
    % RANGE computes the range of the input data X.
    % If X is a vector, it returns the difference between the max and min values.
    % If X is a matrix, it returns a row vector of ranges for each column.
    % If X is a multidimensional array, it operates along the first nonsingleton dimension.
    
    if isempty(X)
        r = []; % Return empty array if input is empty
        return;
    end
    
    % Find the maximum and minimum values along the first nonsingleton dimension
    maxVal = max(X, [], 1);
    minVal = min(X, [], 1);
    
    % Calculate the range
    r = maxVal - minVal;
end