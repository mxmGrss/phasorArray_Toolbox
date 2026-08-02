function [Aph, h1, h2, method] = prepToeplitzInput(Aph, m, method)
%PREPTOEPLITZINPUT  Shared prologue of the Toeplitz assembly kernels.
%
%   [Aph, h1, h2, method] = prepToeplitzInput(Aph, m, method)
%
%   Common to array2TBlocks, array2BToeplitz and sparray2TBlocks, which differ
%   only in how they assemble the final matrix. Unwraps a PhasorArray, resolves
%   m into [h1 h2], and pads/truncates the harmonic dimension to [-h_req,h_req]
%   with h_req = h1+h2 (the widest harmonic their Cauchy product can reach), so
%   harmonic k lands at slice k + h_req + 1.
%
%   m: [] keeps the input order, a scalar is the square case, [h1 h2] is
%   rectangular. method is returned forced to 'cat' for sym/sdpvar inputs,
%   which cell2mat cannot hold.
%
%   See also: array2TBlocks, array2BToeplitz, sparray2TBlocks, isFunny

arguments
    Aph
    m                = []
    method {mustBeMember(method, {'cell2mat', 'cat'})} = 'cell2mat'
end

% --- 1/2. Special-type detection and PhasorArray unwrapping ---------------
% isFunny is checked both on the wrapper and on the unwrapped data: a
% PhasorArray is not itself "funny", its .value may be.
is_special = isFunny(Aph);
if isa(Aph, 'PhasorArray')
    data       = Aph.value;
    is_special = is_special || isFunny(data);
    Aph        = data;
end
if is_special
    method = 'cat';   % cell2mat cannot hold sym/sdpvar objects
end

% --- 3. Resolve the output orders [h1, h2] -------------------------------
n1     = size(Aph, 1);
n2     = size(Aph, 2);
nh_in  = (size(Aph, 3) - 1) / 2;

if isempty(m)
    h1 = nh_in;  h2 = nh_in;
elseif isscalar(m)
    h1 = m;      h2 = m;
else
    h1 = m(1);   h2 = m(2);
end

% --- 4. Pad/truncate the harmonic range to exactly [-h_req, h_req] -------
% h_req = h1+h2 is the widest harmonic a Cauchy product of these orders can
% produce; anything beyond it cannot influence the output blocks.
h_req = h1 + h2;
k_min = -h_req;
k_max =  h_req;

overlap_k_min = max(k_min, -nh_in);
overlap_k_max = min(k_max,  nh_in);

if overlap_k_min <= overlap_k_max
    data_sliced = Aph(:, :, (nh_in + 1 + overlap_k_min) : (nh_in + 1 + overlap_k_max));
    n_lead      = overlap_k_min - k_min;
    n_trail     = k_max - overlap_k_max;
    % cat rather than padarray/indexing so sym/sdpvar survive the padding.
    Aph = cat(3, zeros(n1, n2, n_lead), data_sliced, zeros(n1, n2, n_trail));
else
    % Requested range lies entirely outside the input support.
    Aph = zeros(n1, n2, 2 * h_req + 1);
end
end
