function  barsurf(X,llim,varg)
%BARSURF Create a customizable 3D bar plot.
%
%   BARSURF(X, llim, varg) generates a 3D bar plot of the matrix `X` with
%   options for clipping values, customizing axis labels, scaling, and 
%   adding titles or colorbars.
%
%   Inputs:
%     X          - (matrix) Data to plot as a 3D bar graph.
%     llim       - (scalar) Lower limit for values. Default: 1e-3.
%     Optional name-value parameters:
%                   'title'       - (string) Plot title. Default: ''.
%                   'xticklabel'  - (cell array) X-axis labels. Default: ''.
%                   'yticklabel'  - (cell array) Y-axis labels. Default: ''.
%                   'scale'       - (string) 'log' (default) or 'linear'.
%                   'colorbar'    - (logical) Show colorbar. Default: false.
%                   'EdgeColor'   - (string) Bar edge color. Default: 'black'.
%
%   Outputs:
%     None. Modifies the current figure to display the plot.
%
%   Behavior:
%     - Values in `X` below `llim` are clipped to `llim` for better visualization.
%     - The bar color scales with the `z` values, either logarithmically or linearly.
%     - Custom x/y tick labels and a colorbar can be added for context.
%
%   Example Usage:
%     % Generate a 3D bar plot with linear scaling and a custom title
%     X = rand(5);
%     barsurf(X, 1e-3, struct('title', 'My 3D Bar Plot', 'scale', 'linear'));
%
%     % Add x/y-axis tick labels and enable colorbar
%     barsurf(X, 1e-2, struct('xticklabel', {'A', 'B', 'C', 'D', 'E'}, ...
%                             'yticklabel', {'W', 'X', 'Y', 'Z'}, ...
%                             'colorbar', true));
%
%   See also: bar3, colorbar, gca, title, set.
arguments
    X
    llim = 1e-3
    varg.title=[]
    varg.xticklabel=[]
    varg.yticklabel=[]
    varg.scale {mustBeMember(varg.scale,{'log','linear'})} ='log'
    varg.colorbar=false
    varg.EdgeColor = 'black'
end

if ndims(X)==3
    X=PhasorArray(X);
    X=X.T_tb();
end

if ismatrix(varg.xticklabel) && ~isempty(varg.xticklabel)
    varg.xticklabel=num2cell(varg.xticklabel);
end
if ismatrix(varg.yticklabel) && ~isempty(varg.yticklabel)
    varg.yticklabel=num2cell(varg.yticklabel);
end

X=abs(X);

if isempty(llim)
    llim=min(X,[],"all");
end

% X(X<=llim)=llim;
b=bar3(X);
assignin('base','b',b)
% colormap spring
hchild = get(gca,'Children');
for i = 1:length(hchild)
    ZData = get(hchild(i), 'ZData');
    ZData(ZData<=llim) = llim;
    set(hchild(i), 'ZData', ZData);
end
set(gca,'ZScale',varg.scale)
if varg.colorbar
    colorbar 
    set(gca,'ColorScale',varg.scale)
else
    colorbar off
end

for k = 1:length(b)
    zdata = b(k).ZData;
    switch varg.scale
        case 'log'
            b(k).CData = log10(zdata);
        case 'linear'
            b(k).CData = zdata;
    end
    b(k).FaceColor = 'interp';
    b(k).EdgeColor = varg.EdgeColor;
    set(hggetbehavior(b(k), 'Datacursor'), 'Enable', true);
%     dataTipInteraction('SnapToDataVertex','off')
%     setinteractionhint(b(k), 'DataCursor', true);
end
%     fdgsgh

if ~isempty(varg.title)
    title(varg.title)
end

if ~isempty(varg.xticklabel)
    set(gca,'xtick',1:numel(varg.xticklabel),'xticklabel',[varg.xticklabel{:}])
end

if ~isempty(varg.yticklabel)
    set(gca,'ytick',1:numel(varg.yticklabel),'yticklabel',[varg.yticklabel{:}])
end
a=gca;
a.Interactions = [dataTipInteraction panInteraction];
end

