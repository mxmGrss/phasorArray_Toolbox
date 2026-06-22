function barsurf(o1,thres,hdel,varg)
    % BARSURF Generate a 3D bar surface plot of the phasors of a `PhasorArray`.
    %
    %   BARSURF(o1, thres, hdel, varg) produces a 3D bar surface plot representing
    %   the phasors of `o1`. Phasors below a specified threshold are truncated, except
    %   for the first `hdel` harmonics that meet the condition, which are retained as a margin.
    %
    %   Inputs:
    %     o1     - (PhasorArray) The PhasorArray object to be plotted.
    %     thres  - (double, optional) The relative threshold for truncating phasors.
    %                - Default: `1e-6`.
    %     hdel   - (integer, optional) The number of harmonics to retain as a margin.
    %                - Default: `3`.
    %     varg   - Name-value pair arguments:
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
    %     barsurf(o1);
    %
    %     % Use a custom threshold and linear scaling
    %     barsurf(o1, 1e-6, 3, 'scale', 'linear');
    %
    %   See also: ReduceArray, bar3.
    arguments
        o1
        thres=1e-12;
        hdel=3
        varg.scale {mustBeMember(varg.scale,{'log','linear'})} ='log'
        varg.title  ='Phasor of Matrix'
    end
    if isa(o1.value,'ndsdpvar')
        o1=value(value(o1));
        boolnan=isnan(o1);
        nnz(boolnan);
        if nnz(boolnan)>0
            warning('PhasorArray:double:nanSdpvar', 'Some sdpvar values are NaN; they have been set to 0.')
            o1(boolnan)=0;
            o1=ReduceArray(o1);
        end
    end
    [nx,nz]=size(o1,[1 2]);
    nh=(size(o1,3)-1)/2;
    o1 = o1.neglect(thres,exclude0Phasor=false,reduceMethod="relative");
    [~,refM,hresM]=ReduceArray(o1,reduceMethod="relative",reduceThreshold=thres,exclude0Phasor=false);
    %refM contient en val absolue le plus grand phasor de chaque
    %composante de o1.

    aM = abs(o1.Value);
    maxM = max(aM,[],'all');
    minM = min(aM(aM>0),[],'all'); %minimum non zero value, after the neglect function, every value under was set to 0
    spreadM = maxM/minM;


    %any value lower than minM_signif is considered insignificant
    %minM_signif=min(abs(refM),[],'all')*thres;

    %maxM = max(abs(o1.value),[],'all');
    %maxM = max(abs(o1.value),[],'all');

    %spreadM = maxM/minM_signif;
    logspreadM = log10(spreadM);


    epsM=10^(floor(log10(minM)))*10^-(logspreadM/6);

    hdel=min(nh-hresM,hdel);

    reshM=abs(reshape(o1.value,nx*nz,[]));
    barsurf(reshM(:,((end+1)/2):end).',epsM,"yticklabel",(0:o1.h)','scale',varg.scale)


    %     barsurf(reshM(:,hres+1:end),min(ref,[],'all')*thres,"xticklabel",(0:hres)','scale','log')
    xlabel("States")
    ylabel("Harmonics")
    zlabel('Abs')
    title(varg.title)
    xlim([0 nx*nz+1])

end
