function [Adetph,Adet_t] = PhasorDet(Aph,varg)
    %PHASORDET Computes the phasors of the determinant of a time-varying matrix A(t)
    %
    %   [Adetph, Adet_t] = PhasorDet(Aph, varg) calculates the phasors for the
    %   determinant of a matrix A(t) represented in phasor form, where the input
    %   matrix A(t) is assumed to be periodic in nature. The determinant is computed
    %   for each time slice, and the result is transformed back into phasor form.
    %
    %   Inputs:
    %       - Aph: The 3D array representing the phasors of A(t). The third dimension
    %         corresponds to the harmonic components of A(t).
    %       - varg: Optional arguments, provided as a structure or name-value pairs:
    %           - nT: The number of time steps (default is 1).
    %           - T: The period of the matrix A(t) (default is 1).
    %           - m: The truncation order (default is calculated from the harmonic length).
    %           - plot: Boolean flag to plot results (default is false).
    %           - reduceThreshold: Threshold for reducing the phasors (default is 1e-15).
    %           - reduceMethod: The method to use for reduction (default is 'relative').
    %           - autoTrunc: Flag to automatically truncate based on filtered differences (default is false).
    %
    %   Outputs:
    %       - Adetph: The phasor representation of the determinant of A(t) in 3D array form.
    %       - Adet_t: The time-domain representation of the determinant.
    %
    %   Notes:
    %       - If `Aph` is a `PhasorArray`, it is first converted to its value form.
    %       - The determinant of each time slice of A(t) is computed using the `det` function.
    %       - If `autoTrunc` is set to true, the phasors are truncated based on filtered differences.
    %
    %   Example usage:
    %       [Adetph, Adet_t] = PhasorDet(Aph, 'nT', 10, 'T', 1, 'm', 8, 'plot', true);
    %
    %   See also: PhasorArray, ReduceArray, PhasorArray2time, TimeArray2Phasors

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

nT=varg.nT;
T=varg.T;
m=varg.m;

hA=(size(Aph,3)-1)/2;

if isempty(m)
m=max(ceil(log2(hA+1)+1),8);
end

if hA==0
    Adetph=det(Aph);
    return
end

n=2^m;
t=0:T/n:nT*T-T/n;
if isa(Aph,"ndsdpvar") || isa(Aph,"sdpvar")
    Aph=value(Aph);
end

Aph=ReduceArray(Aph);

At=PhasorArray2time(Aph,T,t,"plot",false);

Adet_t=arrayfun(@(k) det(At(:,:,k)), 1:size(At,3));
Adet_t=permute(Adet_t,[1 3 2]);



Adetphcomp=TimeArray2Phasors(Adet_t,nT);

h=(numel(Adetphcomp)-1)/2;
A_half=squeeze(abs(Adetphcomp(1,1,(h+1):end)));
log10Ph=log10(A_half);
log10Ph(isnan(log10Ph))=mean(log10Ph(end-10:end));
log10Ph(isinf(log10Ph))=mean(log10Ph(end-10:end));
diffPh = log10(abs(diff(A_half,1)));
filt_diff= lowpass(diffPh,0.05);

if ~varg.autoTrunc
Adetph=ReduceArray(Adetphcomp,'reduceMethod',varg.reduceMethod,'reduceThreshold',varg.reduceThreshold);
else
Adetph=ReduceArray(Adetphcomp,find(filt_diff>0,1)+5);
end

if varg.plot
    clf
    tiledlayout("flow")
    nexttile
    h=(numel(Adetph)-1)/2;
    plot(t,squeeze(reshape(Adet_t,[],1,numel(t))))
    hold on
    Ait=PhasorArray2time(Adetph,T,t,"plot",false);
    plot(t,squeeze(reshape(Ait,[],1,numel(t))),'--')
    hold off
    legend('real det','reconstructed det from phasor')
    title('det(A)(t)')
    nexttile
    stem(0:h,squeeze(abs(Adetph(1,1,(h+1):end))))
    set(gca,'YScale','log')
    title('Stem of abs of phasor of det(A)')
    if varg.autoTrunc
    nexttile
    plot(diffPh)
    hold on 
    plot(filt_diff)
    fprintf("recommended truncature order is %d\n",find(filt_diff>0,1))
    title('filtered diff of log10(abs(phasor))')
    end
    
    nexttile
    h=(numel(Adetphcomp)-1)/2;
    stem(0:h,squeeze(abs(Adetphcomp(1,1,(h+1):end))))
    set(gca,'YScale','log')
    title('Stem of abs of phasor of det(A) before trunc')
end

