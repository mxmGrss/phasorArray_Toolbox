function barsurf(pA1,thres,hdel,nvp)
    % BARSURF Generate a 3D bar surface plot of the phasors of a `PhasorArray`.
    %
    %   BARSURF(pA1, thres, hdel, nvp) produces a 3D bar surface plot representing
    %   the phasors of `pA1`. Phasors below a specified threshold are truncated, except
    %   for the first `hdel` harmonics that meet the condition, which are retained as a margin.
    %
    %   Inputs:
    %     pA1     - (PhasorArray) The PhasorArray object to be plotted.
    %     thres  - (double, optional) The relative threshold for truncating phasors.
    %                - Default: `1e-6`.
    %     hdel   - (integer, optional) The number of harmonics to retain as a margin.
    %                - Default: `3`.
    %     nvp   - Name-value pair arguments:
    %                - 'scale' (char): Scale of the plot.
    %                    - 'log' (default): Logarithmic scale.
    %                    - 'linear': Linear scale.
    %
    %   Behavior:
    %     - **Thresholding**: Phasors with absolute values below `thres * maxPhasor`
    %       are truncated, except for the first `hdel` harmonics.
    %     - **Logarithmic Scaling**: By default, the `log` scale is used for better
    %       visualization of magnitude variations.
    %     - **Matrix Reshaping**: The function reshapes the phasors into a suitable
    %       format for `bar3` visualization.
    %
    %   Example Usage:
    %     % Generate a bar plot with default threshold and margin
    %     barsurf(pA1);
    %
    %     % Use a custom threshold and linear scaling
    %     barsurf(pA1, 1e-6, 3, 'scale', 'linear');
    %
    %   See also: ReduceArray, bar3.
    arguments
        pA1
        thres=1e-12;
        hdel=3
        nvp.scale {mustBeMember(nvp.scale,{'log','linear'})} ='log'
        nvp.title  ='Phasor of Matrix'
    end
    if isa(pA1.value,'ndsdpvar')
        pA1=value(value(pA1));
        boolnan=isnan(pA1);
        nnz(boolnan);
        if nnz(boolnan)>0
            warning('PhasorArray:double:nanSdpvar', 'Some sdpvar values are NaN; they have been set to 0.')
            pA1(boolnan)=0;
            pA1=ReduceArray(pA1);
        end
    end
    [nx,nz]=size(pA1,[1 2]);
    nh=(size(pA1,3)-1)/2;
    pA1 = pA1.neglect(thres,"exclude0Phasor", false,"reduceMethod", "relative");
    [~,refM,hresM]=ReduceArray(pA1,"reduceMethod", "relative","reduceThreshold", thres,"exclude0Phasor", false);
    %refM contient en val absolue le plus grand phasor de chaque
    %composante de pA1.

    aM = abs(pA1.Value);
    maxM = max(aM,[],'all');
    minM = min(aM(aM>0),[],'all'); %minimum non zero value, after the neglect function, every value under was set to 0
    spreadM = maxM/minM;


    %any value lower than minM_signif is considered insignificant
    %minM_signif=min(abs(refM),[],'all')*thres;

    %maxM = max(abs(pA1.value),[],'all');
    %maxM = max(abs(pA1.value),[],'all');

    %spreadM = maxM/minM_signif;
    logspreadM = log10(spreadM);


    epsM=10^(floor(log10(minM)))*10^-(logspreadM/6);

    hdel=min(nh-hresM,hdel);

    reshM=abs(reshape(pA1.value,nx*nz,[]));
    barsurf(reshM(:,((end+1)/2):end).',epsM,"yticklabel",(0:pA1.h)','scale',nvp.scale)


    %     barsurf(reshM(:,hres+1:end),min(ref,[],'all')*thres,"xticklabel",(0:hres)','scale','log')
    xlabel("States")
    ylabel("Harmonics")
    zlabel('Abs')
    title(nvp.title)
    xlim([0 nx*nz+1])

end
