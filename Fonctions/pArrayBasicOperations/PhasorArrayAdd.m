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