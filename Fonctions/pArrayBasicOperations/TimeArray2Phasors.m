function [Aph] = TimeArray2Phasors(At,nT,t,varg)
    %TIMEARRAY2PHASORS Convert a time-dependent matrix to its phasor representation.
    %
    %   TIMEARRAY2PHASORS(At, nT, t, varg) converts a 3D array representing a time-dependent
    %   matrix (with time stored along the third dimension) into its phasor representation. 
    %   The phasors are stored along the third dimension of the output array.
    %
    %   Inputs:
    %     At         - (3D array or vector) Input time-domain matrix. If a vector is provided,
    %                  it is assumed to be a scalar signal.
    %     nT         - (scalar) The number of time periods (default: 1).
    %     t          - (vector, optional) Time steps corresponding to the third dimension of `At`.
    %                  Default: `[]`, and the time steps are inferred.
    %     varg       - (struct) Optional parameters:
    %                   'truncIndex'   - (numeric) Index up to which the phasor is truncated. Default: Inf.
    %                   'isReal'       - (logical) If true (default), the positive and negative phasors are 
    %                                     adjusted to satisfy Ak = conj(A-k).
    %                   'timeDim'      - (numeric) The dimension representing time in the 3D array. Default: 3.
    %                   'procedeWithError' - (logical) If true, allows proceeding with errors for non-power-of-2 data.
    %
    %   Outputs:
    %     Aph        - (3D array) The phasor representation of the input time-domain matrix `At`.
    %   
    %   Description:
    %     This function computes the phasors of a time-dependent matrix, assuming the time steps are
    %     uniformly spaced. The Fourier transform is applied along the third dimension of the input
    %     matrix, and the resulting phasors are returned. The positive and negative phasors are adjusted
    %     to respect the condition Ak = conj(A-k) if `isReal` is set to true.
    %     Additionally, it performs validation checks on the input data, such as ensuring the time series
    %     length is a power of 2.
    %
    %   Example Usage:
    %     % Convert time-domain data to phasors with default parameters
    %     Aph = TimeArray2Phasors(At, 1);
    %
    %     % Convert with a specific truncation index and time dimension
    %     Aph = TimeArray2Phasors(At, 1, [], struct('truncIndex', 10, 'timeDim', 2));
    %
    %   See also: fft, ReduceArray, fftshift.
arguments
    At
    nT=1
    t=[]
    varg.truncIndex {mustBeNumeric(varg.truncIndex)} =Inf
    varg.isReal logical = true
    varg.timeDim  {mustBeLessThanOrEqual(varg.timeDim,3),mustBeGreaterThanOrEqual(varg.timeDim,1)} = 3
    varg.procedeWithError = false
end
%  INPUT VALIDATION
if isvector(At) % vector input: reshape to [1 x 1 x n], timeDim is irrelevant
    At=At(:);
    At=permute(At,[2 3 1]);
elseif varg.timeDim~=3 % permute time dimension to dim 3
    if ~ismatrix(At)
        warning('TimeArray2Phasors:timeDimPermutation', 'Non-matrix input with timeDim=%d: permuting dimension %d with dimension 3.', varg.timeDim, varg.timeDim)
    end
    perm         = 1:3;
    perm([varg.timeDim, 3]) = perm([3, varg.timeDim]);
    At           = permute(At, perm);
end

n=size(At,3);

% Validate time vector length against array size
if ~isempty(t) && numel(t) ~= n
    error('PhasorArray:TimeArray2Phasors:timeMismatch', ...
          'numel(t)=%d does not match size(At,timeDim)=%d.', numel(t), n);
end

% Check that samples per period (n/nT) is a power of 2.
% n = nT * 2^m is valid even if n itself is not a power of 2 (e.g. nT=3, n=1536).
nPerPeriod = n / nT;
if ~(nPerPeriod == floor(nPerPeriod) && nPerPeriod > 0 && bitand(uint32(nPerPeriod), uint32(nPerPeriod)-1) == 0)
    if varg.procedeWithError
        warning('TimeArray2Phasors:notPowerOf2', 'Samples per period n/nT=%g is not a power of 2; output accuracy may be reduced.', nPerPeriod)
    else
        error('TimeArray2Phasors:notPowerOf2', 'Samples per period n/nT=%g is not a power of 2; output would be inaccurate. Ensure n = nT * 2^m.', nPerPeriod)
    end
end

FST=fftshift(fft(At,[],3),3)/n;

fshiftn=((-n/2):(n/2-1))/nT;

[~,I1]=find(mod(fshiftn,1)==0);

if mod(numel(I1),2)==0
    Aph=FST(:,:,I1(2:end));
else
    Aph=FST(:,:,I1(1:end));
end

Aph=ReduceArray(Aph,varg.truncIndex);

h=(size(Aph,3)-1)/2;

% Phase correction when t(1) != 0.
% If the signal is sampled on [t1, t1+T-dt] with t1 != 0, each harmonic k
% carries a phase offset exp(-i*k*omega0*t1). T = n*dt/nT.
if ~isempty(t) && abs(t(1)) > eps
    dt_t   = t(2) - t(1);
    T_t    = n * dt_t / nT;
    omega0 = 2*pi / T_t;
    k_vec  = (-h:h);                             % 1 x (2h+1)
    phase  = exp(-1i * k_vec * omega0 * t(1));   % 1 x (2h+1)
    Aph    = Aph .* reshape(phase, 1, 1, []);    % broadcast over dims 1 and 2
end

% Enforce conjugate symmetry Ak = conj(A-k) for real-valued signals (vectorized).
if varg.isReal
    Aph_pos = (Aph(:,:,h+2:end) + conj(flip(Aph(:,:,1:h), 3))) / 2;
    Aph     = cat(3, conj(flip(Aph_pos,3)), ...
                     real(Aph(:,:,h+1)), ...
                     Aph_pos);
end


end

