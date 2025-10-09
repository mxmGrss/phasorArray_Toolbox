function [o2] = PhasorArrayPad2(o1, delta_h)
% PhasorArrayPad Pad a PhasorArray with zeros
% PhasorArrayPad(A, h) pads the PhasorArray A with h zeros on each side of the third dimension
% PhasorArrayPad(A, [n1 n2 ... hn]) pads the PhasorArray A with n1 zeros on the first dimension, n2 zeros on the second dimension, and hn zeros on each side of the nth dimension

dims = size(o1);
numDims = numel(dims);

if numel(delta_h) == 1
    % Pad only the third dimension by default
    padSizes = zeros(1, numDims);
    padSizes(3) = delta_h;
elseif numel(delta_h) == numDims
    % Pad each dimension as specified
    padSizes = delta_h;
else
    error('Invalid input for delta_h in PhasorArrayPad, delta_h must be a scalar or a vector with the same number of elements as dimensions in o1');
end

if isphasor(o1)
    output_phas = 1;
    o1 = o1.value;
else
    output_phas = 0;
end

% Initialize the padded array
switch class(o1)
    case "double"
        o2 = zeros(dims + 2 * padSizes);
    case {"ndsdpvar", "sdpvar"}
        o2 = ndsdpvar(dims + 2 * padSizes, 'full');
    case "sym"
        o2 = sym(zeros(dims + 2 * padSizes));
end

% Determine the indices for the original array in the padded array
indices = arrayfun(@(dim, pad) (pad + 1):(dim + pad), dims, padSizes, 'UniformOutput', false);

% Assign the original array to the appropriate location in the padded array
o2(indices{:}) = o1;

if output_phas
    o2 = PhasorArray(o2, 'reduce', false);
end

end