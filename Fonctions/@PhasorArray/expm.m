function [out] = expm(Aph,nvp)
    %EXPM Compute the matrix exponential of a T-periodic matrix in phasor form.
    %
    %   out = EXPM(Aph, nvp) computes the matrix exponential of a T-periodic
    %   matrix represented as a PhasorArray. The function first converts the
    %   phasor representation into a time-domain matrix, applies the matrix
    %   exponential at each time step, and then transforms it back into phasor form using IFFT.
    %
    %   Inputs:
    %     Aph  - (PhasorArray or numeric) The input periodic matrix in phasor form.
    %     nvp - (Optional) Name-value pair arguments:
    %         'nT'               (integer, default: 1)      - Number of periods for evaluation.
    %         'T'                (double, default: 2*pi)    - Period of the matrix.
    %         'm'                (integer, default: [])     - Resolution parameter for discretization.
    %         'plot'             (logical, default: false)  - Plot the computed result.
    %         'reduceThreshold'  (double, default: 1e-15)   - Threshold for reducing small phasors.
    %         'reduceMethod'     (char, default: 'relative') - Method for phasor reduction.
    %         'autoTrunc'        (logical, default: false)  - Enable automatic truncation.
    %
    %   Output:
    %     out - (PhasorArray) The computed exponential of the input matrix in phasor form.
    %
    %   Example:
    %     Aph = PhasorArray.random(3, 3, 5);
    %     expA = expm(Aph, 'T', 2*pi, 'nT', 1);
    %
    %   See also: LOGM, PhasorArray2time, TimeArray2Phasors, expm.
arguments
Aph
nvp.nT=1
nvp.T=2*pi
nvp.m=[]
nvp.plot=false
nvp.reduceThreshold = 1e-15
nvp.reduceMethod = 'relative'
nvp.autoTrunc = false
end

if isa(Aph,'PhasorArray')
    Aph=Aph.Value;
end
if isa(Aph,"ndsdpvar") || isa(Aph,"sdpvar")
    Aph=value(Aph);
end

nT = nvp.nT;
T  = nvp.T;
m  = nvp.m;

hA=(size(Aph,3)-1)/2;
if isempty(m)
    m=max(ceil(log2(hA+1)+1),8);
end

n  = 2^m;
dt = T/n;

t=0:dt:(nT*T-dt);


if hA==0
    out=PhasorArray(expm(Aph));
    return
end


Aph=ReduceArray(Aph);

At=PhasorArray2time(Aph,T,t,"plot",nvp.plot);

Aex_t=arrayfun(@(k) expm(At(:,:,k)), 1:size(At,3),'UniformOutput',false);
Aex_t=cat(3,Aex_t{:});

out=TimeArray2Phasors(Aex_t,nT);

% The time grid is much finer than the bandwidth of exp(A), so without this
% the result carries n/2 harmonics whatever the input order.
if nvp.autoTrunc
    out=ReduceArray(out,[],'reduceMethod',nvp.reduceMethod, ...
        'reduceThreshold',nvp.reduceThreshold);
end

out = PhasorArray(out);

end

