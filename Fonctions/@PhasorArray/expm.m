function [out] = expm(Aph,varg)
    %EXPM Compute the matrix exponential of a T-periodic matrix in phasor form.
    %
    %   out = EXPM(Aph, varg) computes the matrix exponential of a T-periodic
    %   matrix represented as a PhasorArray. The function first converts the
    %   phasor representation into a time-domain matrix, applies the matrix
    %   exponential at each time step, and then transforms it back into phasor form using IFFT.
    %
    %   Inputs:
    %     Aph  - (PhasorArray or numeric) The input periodic matrix in phasor form.
    %     varg - (Optional) Name-value pair arguments:
    %         'nT'               (integer, default: 1)      - Number of periods for evaluation.
    %         'T'                (double, default: 1)       - Period of the matrix.
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
varg.nT=1
varg.T=1
varg.m=[]
varg.plot=false
varg.reduceThreshold = 1e-15
varg.reduceMethod = 'relative'
varg.autoTrunc = false
end

if isa(Aph,'PhasorArray')
    Aph=Aph.Value;
end
if isa(Aph,"ndsdpvar") || isa(Aph,"sdpvar")
    Aph=value(Aph);
end

nT = varg.nT;
T  = varg.T;
m  = varg.m;

hA=(size(Aph,3)-1)/2;
if isempty(m)
    m=max(ceil(log2(hA+1)+1),8);
end

n  = 2^m;
dt = T/n;

t=0:dt:(nT*T-dt);


if hA==0
    out=expm(Aph);
    return
end


Aph=ReduceArray(Aph);

At=PhasorArray2time(Aph,T,t,"plot",false);

Aex_t=arrayfun(@(k) expm(At(:,:,k)), 1:size(At,3),'UniformOutput',false);
Aex_t=cat(3,Aex_t{:});

Aex_t_ph=TimeArray2Phasors(Aex_t,nT);

out = Aex_t_ph;

end

