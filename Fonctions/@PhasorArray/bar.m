function T = bar(pA, opt)
%BAR Compare harmonic coefficients using grouped or stacked bars.
%   T = bar(A, B, ..., layout="grouped", labels=["A","B"])
%   plots equally sized periodic matrices, aligning different harmonic orders.
%   INPUTS: one or more finite numeric PhasorArray objects. Options follow stem:
%   scale (default "log"), display ("abs", "real", "imag", "both"),
%   side ("oneSided" or "both"), explosed (one tile per matrix entry),
%   parent (figure, axes or tiledlayout), uniformYLim, labels and layout.
%   Complex time signals always use both sides. Values are raw Fourier
%   coefficients, as in stem: positive harmonics are NOT doubled.
%   Signed real/imaginary plots require scale="linear".
%   order=[0 1 3] selects signed orders, overriding side. [] keeps the
%   usual display. Orders are sorted, deduplicated and zero outside storage.
%   OUTPUT: tiledlayout, except an explicit axes parent for a single panel
%   is returned unchanged. Each call groups its own inputs; hold on does
%   not regroup series from earlier calls. No legend is added unless labels
%   are supplied. With explosed=false, matrix entries become separate series.
%   EXAMPLE: bar(A, neglect(A,1e-3), scale="linear", labels=["A","Reduced"])
%   See also stem, pageEnergy, neglect.
arguments (Repeating)
    pA
end
arguments
    opt.layout (1,1) string {mustBeMember(opt.layout,["grouped","stacked"])} = "grouped"
    opt.scale (1,1) string {mustBeMember(opt.scale,["log","linear"])} = "log"
    opt.display (1,1) string {mustBeMember(opt.display,["abs","real","imag","both"])} = "abs"
    opt.side (1,1) string {mustBeMember(opt.side,["oneSided","both"])} = "oneSided"
    opt.explosed (1,1) logical = true
    opt.parent = []
    opt.uniformYLim (1,1) logical = false
    opt.labels (1,:) string = strings(1,0)
    opt.order {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger} = []
end
if isempty(pA) || ~all(cellfun(@(a) isa(a,'PhasorArray'),pA))
    error('PhasorArray:bar:input','Provide one or more PhasorArray objects.');
end
nr=size(pA{1},1); nc=size(pA{1},2); n=numel(pA);
if nr*nc==0 || ~all(cellfun(@(a) size(a,1)==nr && size(a,2)==nc,pA))
    error('PhasorArray:bar:dimensions','Inputs must have identical, nonempty matrix dimensions.');
end
if ~isempty(opt.labels) && numel(opt.labels)~=n
    error('PhasorArray:bar:labels','Provide one label per PhasorArray.');
end
if opt.scale=="log" && opt.display~="abs"
    error('PhasorArray:bar:signedLog','Use scale="linear" for signed real/imaginary coefficients.');
end
H=max(cellfun(@(a) a.h,pA)); orders=-H:H;
values=zeros(nr,nc,2*H+1,n);
for k=1:n
    v=pvalue(pA{k});
    if ~isnumeric(v) || any(~isfinite(v(:)))
        error('PhasorArray:bar:payload','Plotting requires finite numeric coefficients.');
    end
    values(:,:,H+1+(-pA{k}.h:pA{k}.h),k)=v;
end
if ~isempty(opt.order)
    validateattributes(opt.order,{'numeric'},{'vector'});
    orders=unique(double(opt.order(:).'));
    selected=zeros(nr,nc,numel(orders),n);
    valid=abs(orders)<=H;
    selected(:,:,valid,:)=values(:,:,H+1+orders(valid),:);
    values=selected;
elseif opt.side=="oneSided" && all(cellfun(@(a) isreal(a),pA))
    orders=0:H; values=values(:,:,H+1:end,:);
end
parts=opt.display;
if parts=="both", parts=["real","imag"]; end
if opt.explosed, rows=nr; cols=nc; else, rows=1; cols=1; end
panels=rows*cols;
parent=opt.parent;
if isempty(parent), parent=gcf; end
if isgraphics(parent,'axes') && panels==1
    T=parent; axesList=parent;
else
    if isa(parent,'matlab.graphics.layout.TiledChartLayout')
        T=tiledlayout(parent,rows,cols,'TileSpacing','compact','Padding','compact');
    elseif isgraphics(parent,'figure') || isgraphics(parent,'axes')
        T=manageTiledLayout(parent,rows,cols,"plotBarPhasor");
    else
        error('PhasorArray:bar:parent','Parent must be a figure, axes or tiledlayout.');
    end
    axesList=gobjects(1,panels);
    for panel=1:panels, axesList(panel)=nexttile(T,panel); end
end
for panel=1:panels
    if opt.explosed
        [j,i]=ind2sub([nc,nr],panel); entries=sub2ind([nr,nc],i,j);
    else
        entries=1:nr*nc;
    end
    Y=[]; names=strings(1,0);
    for k=1:n
        v=reshape(values(:,:,:,k),nr*nc,[]);
        for entry=entries
            [i,j]=ind2sub([nr,nc],entry);
            for part=parts
                switch part
                    case "abs", y=abs(v(entry,:));
                    case "real", y=real(v(entry,:));
                    case "imag", y=imag(v(entry,:));
                end
                Y(:,end+1)=y.'; %#ok<AGROW>
                label="A"+k;
                if ~isempty(opt.labels), label=opt.labels(k); end
                if ~opt.explosed && nr*nc>1, label=label+sprintf(' (%d,%d)',i,j); end
                if numel(parts)>1, label=label+" "+part; end
                names(end+1)=label; %#ok<AGROW>
            end
        end
    end
    % MATLAB treats a single-row Y as a vector, not multiple series. An
    % invisible NaN group keeps the DC-only case grouped by input object.
    x=orders;
    if numel(x)==1, x=[x,x+1]; Y=[Y;nan(1,size(Y,2))]; end
    handles=bar(axesList(panel),x,Y,char(opt.layout));
    set(axesList(panel),'YScale',char(opt.scale));
    xlim(axesList(panel),[min(orders)-0.6,max(orders)+0.6]);
    xticks(axesList(panel),orders); grid(axesList(panel),'on');
    xlabel(axesList(panel),'Harmonic order'); ylabel(axesList(panel),opt.display);
    if opt.explosed && nr*nc>1, title(axesList(panel),sprintf('(%d,%d)',i,j)); end
    for k=1:numel(handles), handles(k).DisplayName=names(k); end
    if panel==1 && ~isempty(opt.labels), legend(axesList(panel),'show'); end
end
if panels>1
    linkaxes(axesList,'x');
    if opt.uniformYLim, linkaxes(axesList,'y'); end
end
end
