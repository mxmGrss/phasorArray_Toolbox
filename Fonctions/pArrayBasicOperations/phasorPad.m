function Apadded = phasorPad(A, padsize, padval, direction)
%PHASORPAD Pad a 3D array (PhasorArray compatible replacement for padarray)
%
%   Apadded = phasorPad(A, padsize) pads array A with zeros.
%   Apadded = phasorPad(A, padsize, padval) pads with the specified value.
%   Apadded = phasorPad(A, padsize, padval, direction) specifies padding direction.
%
%   This function provides a MATLAB-native alternative to padarray() from the
%   Signal Processing Toolbox, ensuring compatibility with:
%       - Numeric arrays (double, single)
%       - Symbolic arrays (sym)
%       - YALMIP decision variables (sdpvar, ndsdpvar)
%       - PhasorArray objects
%
%   SYNTAX:
%       Apadded = phasorPad(A, padsize)
%       Apadded = phasorPad(A, padsize, padval)
%       Apadded = phasorPad(A, padsize, padval, direction)
%
%   INPUTS:
%       A        - Input array (numeric, sym, sdpvar, or PhasorArray)
%       padsize  - Padding size specification:
%                  • Scalar: pad all dimensions equally
%                  • [P]: pad all dimensions by P
%                  • [P1 P2 P3]: pad each dimension by [rows, cols, pages]
%       padval   - (Optional) Padding value (default: 0)
%       direction - (Optional) Padding direction:
%                  • 'both' (default): pad on both sides
%                  • 'pre': pad before (left/top/front)
%                  • 'post': pad after (right/bottom/back)
%
%   OUTPUTS:
%       Apadded  - Padded array (same type as input)
%
%   EXAMPLES:
%       % Pad numeric array
%       A = randn(2, 3, 5);
%       Apad = phasorPad(A, [1 2 3]);
%       % Result: size(Apad) = [4, 7, 11]  (both sides)
%
%       % Pad symbolic array
%       Asym = sym('a', [2, 2, 3]);
%       Apad = phasorPad(Asym, [0 0 2], 0, 'post');
%       % Result: size(Apad) = [2, 2, 7]  (post only)
%
%       % Pad PhasorArray (harmonics dimension)
%       A = PhasorArray(randn(3, 3, 7));
%       Apad = phasorPad(A.Value, [0 0 2]);
%       % Result: size(Apad) = [3, 3, 11]
%
%   COMPATIBILITY:
%       This function is fully compatible with:
%       - MATLAB R2021b+ (base MATLAB, no toolboxes required)
%       - Symbolic Math Toolbox (sym arrays)
%       - YALMIP (sdpvar, ndsdpvar)
%
%   See also: cat, zeros, padarray (Signal Processing Toolbox)
%
%   PhasorArray Toolbox
%   Author: Maxime GROSSO
%   Date: March 2026

arguments
    A
    padsize (1,:) {mustBeInteger, mustBeNonnegative}
    padval = 0
    direction {mustBeMember(direction, {'both', 'pre', 'post'})} = 'both'
end

% Extract value if input is PhasorArray
if isa(A, 'PhasorArray')
    A = A.Value;
end

% Get array dimensions
sz = size(A);
ndim = length(sz);

% Parse padsize input
if isscalar(padsize)
    % Scalar: pad all dimensions equally
    padsize = repmat(padsize, 1, ndim);
elseif length(padsize) < ndim
    % Extend padsize with zeros for missing dimensions
    padsize = [padsize, zeros(1, ndim - length(padsize))];
elseif length(padsize) > ndim
    % Truncate or error
    error('phasorPad:InvalidPadsize', ...
          'padsize has %d elements but array has %d dimensions.', ...
          length(padsize), ndim);
end

% Determine padding for each side
switch direction
    case 'both'
        pad_pre = padsize;
        pad_post = padsize;
    case 'pre'
        pad_pre = padsize;
        pad_post = zeros(1, ndim);
    case 'post'
        pad_pre = zeros(1, ndim);
        pad_post = padsize;
end

% Create padding template (zeros or specified value)
% Must handle sym and sdpvar types correctly
if isa(A, 'sym')
    % Symbolic: create sym zeros - use manual construction
    isSymbolic = true;
    isYalmip = false;
elseif isa(A, 'sdpvar') || isa(A, 'ndsdpvar')
    % YALMIP: create sdpvar zeros using arithmetic
    isSymbolic = false;
    isYalmip = true;
else
    % Numeric: standard zeros
    isSymbolic = false;
    isYalmip = false;
end

% Helper function to create zeros of appropriate type
createZeros = @(sz_vec) createZerosHelper(A, sz_vec, isSymbolic, isYalmip);

% Apply padding dimension by dimension
Apadded = A;

% Dimension 1 (rows)
if pad_pre(1) > 0 || pad_post(1) > 0
    sz_current = size(Apadded);
    if pad_pre(1) > 0
        sz_pad = [pad_pre(1), sz_current(2:end)];
        pad_top = createZeros(sz_pad);
        if padval ~= 0
            pad_top = pad_top + padval;
        end
        Apadded = cat(1, pad_top, Apadded);
    end
    if pad_post(1) > 0
        sz_current = size(Apadded);
        sz_pad = [pad_post(1), sz_current(2:end)];
        pad_bottom = createZeros(sz_pad);
        if padval ~= 0
            pad_bottom = pad_bottom + padval;
        end
        Apadded = cat(1, Apadded, pad_bottom);
    end
end

% Dimension 2 (columns)
if length(padsize) >= 2 && (pad_pre(2) > 0 || pad_post(2) > 0)
    sz_current = size(Apadded);
    if pad_pre(2) > 0
        sz_pad = [sz_current(1), pad_pre(2)];
        if length(sz_current) > 2
            sz_pad = [sz_pad, sz_current(3:end)];
        end
        pad_left = createZeros(sz_pad);
        if padval ~= 0
            pad_left = pad_left + padval;
        end
        Apadded = cat(2, pad_left, Apadded);
    end
    if pad_post(2) > 0
        sz_current = size(Apadded);
        sz_pad = [sz_current(1), pad_post(2)];
        if length(sz_current) > 2
            sz_pad = [sz_pad, sz_current(3:end)];
        end
        pad_right = createZeros(sz_pad);
        if padval ~= 0
            pad_right = pad_right + padval;
        end
        Apadded = cat(2, Apadded, pad_right);
    end
end

% Dimension 3 (pages/harmonics)
if length(padsize) >= 3 && (pad_pre(3) > 0 || pad_post(3) > 0)
    sz_current = size(Apadded);
    if pad_pre(3) > 0
        pad_front = createZeros([sz_current(1), sz_current(2), pad_pre(3)]);
        if padval ~= 0
            pad_front = pad_front + padval;
        end
        Apadded = cat(3, pad_front, Apadded);
    end
    if pad_post(3) > 0
        sz_current = size(Apadded);
        pad_back = createZeros([sz_current(1), sz_current(2), pad_post(3)]);
        if padval ~= 0
            pad_back = pad_back + padval;
        end
        Apadded = cat(3, Apadded, pad_back);
    end
end

end

% Helper function to create zeros with proper type
function Z = createZerosHelper(A, sz_vec, isSymbolic, isYalmip)
    if isSymbolic
        % Symbolic: create using sym(zeros(...))
        Z = sym(zeros(sz_vec));
    elseif isYalmip
        % YALMIP: multiply a scalar element by zero and replicate
        Z = 0 * A(1,1,1) * ones(sz_vec);
    else
        % Numeric: standard zeros
        Z = zeros(sz_vec);
    end
end
