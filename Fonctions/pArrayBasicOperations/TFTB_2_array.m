function [Xph] = TFTB_2_array(colX,n1,n2)
%TFTB_2_ARRAY Convert a column vector of stacked TFs (Fourier Transforms) into a 3D array.
%
%   TFTB_2_ARRAY(colX, n1, n2) takes a column vector `colX` composed of stacked 
%   Fourier Transforms (TF) of order `h` of signals, and converts it into a 3D array 
%   representing the original size of the signal `x`.
%
%   The input `colX` is assumed to be a column vector of length `n1 * n2 * (2h + 1)`, 
%   where `n1` and `n2` are the first and second dimensions of `x`, respectively. The 
%   number of harmonics `h` is inferred from the length of `colX` and `n1, n2`.
%
%   Inputs:
%     colX    - (column vector) A column vector of length `n1 * n2 * (2h + 1)`,
%               composed of stacked Fourier Transforms of signals.
%     n1      - (scalar) The first dimension of the original signal `x`. Default: 1.
%     n2      - (scalar, optional) The second dimension of the original signal `x`. 
%               Default: 1. If omitted, the signal is considered to be a vector.
%
%   Outputs:
%     Xph     - (3D array) The original signal in a 3D array of size `n1 x n2 x (2h + 1)`.
%               The third dimension corresponds to the harmonics of the Fourier Transform.
%
%   Description:
%     This function converts a column vector `colX` of stacked Fourier Transforms of 
%     order `h` into a 3D array. The length of `colX` is expected to be `n1 * n2 * (2h + 1)`,
%     where `h` is deduced from the length of `colX` and the dimensions `n1` and `n2`.
%     If `n2` is omitted, it is assumed to be equal to 1, and the signal is considered a vector. 
%     If both `n1` and `n2` are omitted, the signal is considered to be scalar with `n1 = n2 = 1`.
%     The function reshapes the input vector into a 3D array with the given dimensions, 
%     and the third dimension contains the harmonics of the signal.
%
%   Example Usage:
%     % Convert a column vector of stacked Fourier Transforms into a 3D array
%     Xph = TFTB_2_array(colX, 3, 3);
%
%     % Convert a column vector of stacked Fourier Transforms for a scalar signal
%     Xph = TFTB_2_array(colX);
%
%   See also: reshape, permute, error.
arguments
    colX
    n1=1
    n2=1
end

if size(colX,2)>1
    if n2 == 1
        n2 = size(colX,2);
    elseif n2 ~= size(colX,2)
        error("n2 should be equal to second dim of colX when colX is matricial")
    end

    for iterCol_i = 1:size(colX,2)
        XCol_i = colX(:,iterCol_i);
        % outX(:,iterCol_i,:) = TFTB_2_array(XCol_i,n1,1);
        outX{iterCol_i} = TFTB_2_array(XCol_i,n1,1);
    end
    Xph = cat(2,outX{:});
    return
end



if isrow(colX)
    colX=colX.';
end

m=numel(colX);
dHp1=m/n1/n2;
h=(dHp1-1)/2;

if h~=round(h)
    error('Wrong dimension, non integer number of harmonics, pls check arguments')
end

dX=reshape(full(colX),dHp1,n1,n2);
Xph=permute(dX,[2,3,1]);
end

