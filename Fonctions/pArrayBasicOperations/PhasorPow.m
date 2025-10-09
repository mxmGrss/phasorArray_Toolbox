function [Apph] = PhasorPow(Aph,p,varg)
    % PHASORPOW Computes the matrix power \( A(t)^p \) for a periodic matrix \( A(t) \).
    %
    %   Apph = PhasorPow(Aph, p) computes the matrix power \( A(t)^p \) at each time instant,
    %   reconstructs its phasor representation, and returns it as a PhasorArray.
    %
    %   Inputs:
    %   - Aph: A PhasorArray (or a 3D array representing phasors) of a periodic matrix \( A(t) \).
    %   - p: The power exponent to apply (default: 1).
    %   - varg (Name-Value Arguments):
    %       - nT (default: 1): Number of periods considered in the computation.
    %       - T (default: 1): The fundamental period of \( A(t) \).
    %       - m (default: []): Controls the discretization size; if empty, an automatic choice is made.
    %       - plot (default: false): If true, compares time-domain results before and after phasor conversion.
    %
    %   Output:
    %   - Apph: The phasor representation of \( A(t)^p \), returned as a PhasorArray.
    %
    %   Behavior:
    %   - Converts \( A(t) \) from phasor representation to time domain.
    %   - Computes \( A(t)^p \) pointwise in the time domain.
    %   - Converts the result back into phasor form.
    %   - If `plot = true`, it visualizes both the computed power and its phasor approximation.
    %
    %   Notes:
    %   - This corresponds to the **true matrix power** \( A(t)^p \) in the time domain.
    %   - For non-integer \( p \), matrix functions like `logm` and `expm` may be involved.
    %   - This function does **not** apply exponentiation directly to phasors.
    %
    %   See also: power, mpower, logm, expm, PhasorArray2time, TimeArray2Phasors

arguments
Aph
p=1
varg.nT=1
varg.T=1
varg.m=[]
varg.plot=false
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

n=2^m;
t=0:T/n:nT*T-T/n;
Aph=ReduceArray(Aph);
At=PhasorArray2time(Aph,T,t,"plot",false);

% t1=tic();
% Atp=arrayfun(@(ii) At(:,:,ii)^p,1:n,'UniformOutput',false);
% Atp=cat(3,Atp{:});
% 
% toc(t1)
% t1=tic();

Atp=At*0;
for ii=1:n
    Atp(:,:,ii)=At(:,:,ii)^p;
end

% toc(t1)

Apph=TimeArray2Phasors(Atp,nT,t,isReal=false);

if varg.plot
    plot(t,squeeze(reshape(Atp,[],1,numel(t))))
    hold on
    Aipt=PhasorArray2time(Apph,T,t,"plot",false);
    plot(t,squeeze(reshape(Aipt,[],1,numel(t))),'--')
    hold off
end

