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
% dt=t(2)-t(1);
% T=(t(end)+dt)/nT;


%  VALIDATION DE L'ENTREE
if isvector(At) %case of vector input, dont need to have a 3D structure, dont need to specify timeDim
    At=At(:); % At is put as a column
    At=permute(At,[2 3 1]);
elseif varg.timeDim~=3 %otherwise we check if timeDim is specified
    if ismatrix(At) %Case of matrix is provided
        %all good
    else
        %send a warning, does not respect convention
        warning('specified timeDim to value ~=3 with non matrix input(3D Array), permuting %d dim with 3rd dim',varg.timeDim)
    end
    switch varg.timeDim %perform permutation
        case 1
            At = permute(At,[3 2 1]);
        case 2
            At = permute(At,[1 3 2]);
    end
end

n=size(At,3);

testPowerOf2 = log(n)/log(2);
if abs(round(testPowerOf2)-testPowerOf2)>1e-3
    if varg.procedeWithError
        warning('it looks like your time series is not of length of the form 2^m, ouput will be unacurrate')
    else
        error('it looks like your time series is not of length of the form 2^m (length = %d), ouput would be unacurrate, recheck your data or pass the argument procedeWithError to true',n)
    end
end

FST=fftshift(fft(At,[],3),3)/n;

fshiftn=((-n/2):(n/2-1))/nT;
% fshift=fshiftn/T;

[~,I1]=find(mod(fshiftn,1)==0);

if mod(numel(I1),2)==0
    Aph=FST(:,:,I1(2:end));
else
    Aph=FST(:,:,I1(1:end));
end

Aph=ReduceArray(Aph,varg.truncIndex);

h=(size(Aph,3)-1)/2;
if varg.isReal
    for kk=0:h
        Aph(:,:,h+1+kk)=(Aph(:,:,h+1+kk)+conj(Aph(:,:,h+1-kk)))/2;
        Aph(:,:,h+1-kk)=conj(Aph(:,:,h+1+kk));
    end
end


end

