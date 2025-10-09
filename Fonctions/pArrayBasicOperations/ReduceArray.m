function [Aphr,ref,htrunc] = ReduceArray(Aph,htrunc,varg)
%REDUCEARRAY Reduce the PhasorArray according to various methods
%   [Aphr, ref, htrunc] = REDUCEARRAY(Aph) reduces the PhasorArray Aph by truncating excessive phasors.
%   [Aphr, ref, htrunc] = REDUCEARRAY(Aph, htrunc) truncates the PhasorArray to the specified harmonic order htrunc.
%   [Aphr, ref, htrunc] = REDUCEARRAY(Aph, htrunc, 'reduceMethod', method, 'reduceThreshold', threshold, 'exclude0Phasor', exclude, 'hardThresholdPhasors', hardThreshold)
%   reduces the PhasorArray according to the specified method and threshold.
%
%   Input Arguments:
%   Aph - The PhasorArray object or 3D array of phasors to be reduced.
%   htrunc - (Optional) The harmonic order to truncate the PhasorArray to.
%   varg - (Optional) Name-value pair arguments:
%       'reduceMethod' - The method to use for reduction ('absolute' or 'relative'). Default is 'absolute'.
%       'reduceThreshold' - The threshold value for reduction. Default is 0.
%       'exclude0Phasor' - Logical flag to exclude the 0th phasor from consideration. Default is false.
%       'hardThresholdPhasors' - Logical flag to apply a hard threshold to phasors. Default is false.
%
%   Output Arguments:
%   Aphr - The reduced PhasorArray.
%   ref - The reference phasor used for relative reduction.
%   htrunc - The harmonic order to which the PhasorArray was truncated.
%
%   Example:
%   [Aphr, ref, htrunc] = ReduceArray(Aph);
%   [Aphr, ref, htrunc] = ReduceArray(Aph, 5);
%   [Aphr, ref, htrunc] = ReduceArray(Aph, [], 'reduceMethod', 'relative', 'reduceThreshold', 1e-10, 'exclude0Phasor', true);
%
%   See also: PhasorArray
%
arguments
    Aph
    htrunc=[]
    varg.reduceMethod char {mustBeMember(varg.reduceMethod,{'absolute','relative'})} = 'absolute'
    varg.reduceThreshold {mustBeNumeric,mustBeReal} = 0
    varg.exclude0Phasor (1,1) logical = false
    varg.hardThresholdPhasors = false
end


if isa(Aph,'PhasorArray')
    Aph=Aph.Value;
end

nhA=(size(Aph,3)-1)/2;

if ~isempty(htrunc)
    htrunc=min(htrunc,nhA);
    Aphr=Aph(:,:,(nhA+1)+(-htrunc:htrunc));
    return
else

    %maximal ref phasor for each coordinate
    if varg.exclude0Phasor
        ref=max(abs(Aph(:,:,[1:nhA , (nhA+2):(2*nhA+1)])),[],3); %maximum harmonic on each coeef, excepting the phasor 0.
    else
        ref=max(abs(Aph),[],3); %maximum harmonic on each coeef.
    end

    kk=nhA;
    switch varg.reduceMethod
        case 'absolute'
            toto = abs(Aph(:,:,nhA+1-kk)) + abs(Aph(:,:,nhA+1+kk));
            while isempty(find(toto>varg.reduceThreshold,1)) && kk>0
                kk=kk-1;
                toto = abs(Aph(:,:,nhA+1-kk)) + abs(Aph(:,:,nhA+1+kk));
            end
            Aphr=Aph(:,:,nhA+1+(-kk:kk));
            htrunc=kk;


            if varg.hardThresholdPhasors
                Aphr(abs(Aphr) < varg.reduceThreshold) = 0;
            end


        case 'relative'
            toto = abs(Aph(:,:,nhA+1-kk)) + abs(Aph(:,:,nhA+1+kk));
            toto = toto./2./ref;
            while isempty(find(toto>varg.reduceThreshold,1)) && kk>0
                kk=kk-1;
                toto = abs(Aph(:,:,nhA+1-kk)) + abs(Aph(:,:,nhA+1+kk));
                toto = toto./2./ref;
            end
            Aphr=Aph(:,:,nhA+1+(-kk:kk));
            htrunc=kk;

            if varg.hardThresholdPhasors
                Aphr(abs(Aphr)./ref < varg.reduceThreshold) = 0;
            end
    end
end
end