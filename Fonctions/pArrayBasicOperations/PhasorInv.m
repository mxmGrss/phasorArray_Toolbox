function [Ainvph,At,norm_err,norm_ref] = PhasorInv(Aph,varg)
%PHASORINV Compute the phasor representation of A(t)^-1 via pointwise inversion.
%
%   This function computes the **pointwise inverse** of the time-domain matrix A(t)
%   and reconstructs its phasor representation.
%
%   Steps:
%   1. Convert phasor representation Aph into the time-domain signal A(t) via **IFFT**.
%   2. Compute A⁻¹(t) **pointwise** for each time sample.
%      - If A(t) is square, uses **pageinv()** for direct inversion.
%      - If A(t) is non-square, uses **pseudo-inversion** via pinv().
%   3. Reconstruct phasors by applying an **FFT** to A⁻¹(t).
%   4. Optionally truncates phasors based on reduction methods.
%
%   Syntax:
%   [Ainvph, At, norm_err, norm_ref] = PHASORINV(Aph)
%       Computes the phasors of A⁻¹(t) using default settings.
%
%   [Ainvph, At, norm_err, norm_ref] = PHASORINV(Aph, 'nT', nT, 'T', T, 'm', m, 
%                                                      'plot', plotFlag, 'autoTrunc', autoTrunc)
%       Computes A⁻¹(t) with additional control over truncation, thresholding, and plotting.
%
%   Input Arguments:
%   - Aph (PhasorArray | numeric array) : The phasor representation of A(t).
%
%   Name-Value Pair Arguments:
%   - 'nT' (integer, optional) : Number of periods used in the time-domain evaluation. Default: 1.
%   - 'T' (double, optional) : The period used for simulation. Default: 1.
%   - 'm' (integer, optional) : 
%       - Power of two controlling time-domain discretization.
%       - Can be set to [] for automatic selection based on the number of phasors.
%   - 'plot' (logical, optional) : If true, plots A⁻¹(t) after computation. Default: false.
%   - 'autoTrunc' (logical, optional) : 
%       - true : Uses the derivative of phasors to **automatically detect** 
%         the significant number of phasors.
%       - false (default) : Uses a fixed threshold-based reduction method.
%
%   - If 'autoTrunc' is false, the following options apply:
%       - 'reduceThreshold' (double, optional) : The threshold for reducing phasors. Default: 4e-15.
%       - 'reduceMethod' (char, optional) : Reduction strategy.
%           - 'absolute' : Remove phasors with magnitude < reduceThreshold.
%           - 'relative' (default) : Remove phasors with magnitude < max(magnitude) * reduceThreshold.
%   - 'verbose' (logical, optional) : If true, prints debug information. Default: false.
%
%   Output Arguments:
%   - Ainvph (PhasorArray) : The phasor representation of A⁻¹(t).
%   - At (array) : The time-domain realization of A(t).
%   - norm_err (double) : Reconstruction error ||Ainv_ph - Ainv_t||_F.
%   - norm_ref (double) : Reference norm ||Ainv_t||_F.
%
%   Example:
%   % Compute the phasors of A⁻¹(t) using default settings
%   [Ainvph, At] = PhasorInv(A);
%
%   % Compute A⁻¹(t) over 2 periods with auto truncation
%   [Ainvph, At] = PhasorInv(A, 'nT', 2, 'autoTrunc', true);
%
%   % Compute A⁻¹(t) with manual truncation and thresholding
%   [Ainvph, At] = PhasorInv(A, 'reduceThreshold', 1e-15, 'reduceMethod', 'absolute');
%
%   See also: DET, REDUCE, FFT, IFFT.

arguments
    Aph
    varg.nT=1
    varg.T=1
    varg.m=[]
    varg.plot=false
    varg.reduceThreshold = 4e-15
    varg.reduceMethod = 'relative'
    varg.autoTrunc = false
    varg.verbose = false
    varg.evalInv = false
end

if isa(Aph,'PhasorArray')
    Aph=Aph.Value;
end

nT=varg.nT;
T=varg.T;
m=varg.m;

hA=(size(Aph,3)-1)/2;

m_nyquist  = nextpow2(hA*2);
m_reco_min = max(nextpow2(hA*4),8);

if isempty(m)
    m=m_reco_min;
end
if varg.verbose
    disph("Minimal m for Nyquist is mN = ",m_nyquist, ", recommanded m is mN + 1 =  ", m_reco_min,", used m is ", m)
end

if hA==0
    Ainvph=Aph^-1;
    return
end

n=2^m;
t=0:T/n:nT*T-T/n;
if isa(Aph,"ndsdpvar") || isa(Aph,"sdpvar")
    Aph=value(Aph);
end

Aph=ReduceArray(Aph);

At=PhasorArray2time(Aph,T,t,"plot",false);
try 
    Ainvt=pageinv(At);
    catch e1
        warning(e1.message)
        warning("pageinv failed, using pagepinv")
    try
        Ainvt=pagepinv(At);
    catch e2
        warning(e2.message)
        warning("pagepinv failed, using pinv")
        for ii = 1:size(At,3)
            Ainvt(:,:,ii)=pinv(At(:,:,ii));
        end
    end
end
Ainvph_compu=TimeArray2Phasors(Ainvt,nT);
if varg.autoTrunc
    Ainvph=ReduceArray(Ainvph_compu,'reduceMethod',varg.reduceMethod,'reduceThreshold',varg.reduceThreshold);
else
    % Ainvph=ReduceArray(Ainvph_compu,find(filt_diff>0,1)+5);
    Ainvph=Ainvph_compu;
end

h=(size(Ainvph,3)-1)/2;
log10Ph=log10(abs(Ainvph(:,:,(h+1):end)));
toto=numel(find(isnan(log10Ph)));
if toto>0
    log10Ph(isnan(log10Ph))=repmat(mean(log10Ph(:,:,end-10:end),'all'),[1 1 toto]);
end
toto=numel(find(isinf(log10Ph)));
if toto>0
    log10Ph(isinf(log10Ph))=repmat(mean(log10Ph(:,:,end-10:end),'all'),[1 1 toto]);
end
% diffPh = diff(log10Ph,1,3);
% filt_diff= lowpass(diffPh,0.05);

if varg.verbose || varg.plot || varg.evalInv || nargout>2
    Ait=PhasorArray2time(Ainvph,T,t,"plot",false);
    Err_recons=Ait-Ainvt;
    norm_err=norm(Err_recons,"fro");
    norm_ref=norm(Ainvt,"fro");
    if varg.verbose
    disph("global error energy is ",norm_err/norm_ref)
    end
end

if varg.plot
    TL=tiledlayout("flow");
    TTL1=tiledlayout(TL,size(log10Ph,1),size(log10Ph,2));
    TTL1.Layout.Tile=1;
    for ii = 1:size(log10Ph,1)
        for jj = 1:size(log10Ph,2)
            nexttile(TTL1)
            plot(t,real(squeeze(Ainvt(ii,jj,:))))
            hold on
            plot(t,real(squeeze(reshape(Ait(ii,jj,:),[],1,numel(t)))),'--')
            hold off
            grid on
            grid minor
            title(sprintf("error energy is %d",norm(squeeze(Err_recons(ii,jj,:)))))
        end
    end
    % plot(t,squeeze(reshape(Ainvt,[],1,numel(t))))
    % hold on
    % Ait=PhasorArray2time(Ainvph,T,t,"plot",false);
    % plot(t,squeeze(reshape(Ait,[],1,numel(t))),'--')
    % hold off
    % grid on
    % grid minor
    TTL=tiledlayout(TL,size(log10Ph,1),size(log10Ph,2));
    TTL.Layout.Tile=2;
    for ii = 1:size(log10Ph,1)
        for jj = 1:size(log10Ph,2)
            nexttile(TTL)
            stem(0:h,squeeze(abs(Ainvph(ii,jj,(h+1):end))))
            set(gca,'YScale','log')
        end
    end



end

