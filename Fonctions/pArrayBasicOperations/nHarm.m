function n = nHarm(x)
%NHARM Number of positive harmonics of a phasor array stored along dim 3.
%   For a raw array X of size m-by-n-by-(2h+1), returns h = (size(X,3)-1)/2.
%   Operates on the raw 3D array (use the .h property for a PhasorArray).
arguments
    x
end
n = (size(x, 3) - 1) / 2;
end
