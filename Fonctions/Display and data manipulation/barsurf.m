function  barsurf(X,llim,nvp)
%BARSURF Create a customizable 3D bar plot.
%
%   BARSURF(X, llim, nvp) generates a 3D bar plot of the matrix `X` with
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
    nvp.title=[]
    nvp.xticklabel=[]
    nvp.yticklabel=[]
    nvp.scale {mustBeMember(nvp.scale,{'log','linear'})} ='log'
    nvp.colorbar=false
    nvp.EdgeColor = 'black'
end

if ndims(X)==3
    X=PhasorArray(X);
    X=X.T_tb();
end

if ismatrix(nvp.xticklabel) && ~isempty(nvp.xticklabel)
    nvp.xticklabel=num2cell(nvp.xticklabel);
end
if ismatrix(nvp.yticklabel) && ~isempty(nvp.yticklabel)
    nvp.yticklabel=num2cell(nvp.yticklabel);
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
set(gca,'ZScale',nvp.scale)
if nvp.colorbar
    colorbar 
    set(gca,'ColorScale',nvp.scale)
else
    colorbar off
end

for k = 1:length(b)
    zdata = b(k).ZData;
    switch nvp.scale
        case 'log'
            b(k).CData = log10(zdata);
        case 'linear'
            b(k).CData = zdata;
    end
    b(k).FaceColor = 'interp';
    b(k).EdgeColor = nvp.EdgeColor;
    set(hggetbehavior(b(k), 'Datacursor'), 'Enable', true);
%     dataTipInteraction('SnapToDataVertex','off')
%     setinteractionhint(b(k), 'DataCursor', true);
end
%     fdgsgh

if ~isempty(nvp.title)
    title(nvp.title)
end

if ~isempty(nvp.xticklabel)
    set(gca,'xtick',1:numel(nvp.xticklabel),'xticklabel',[nvp.xticklabel{:}])
end

if ~isempty(nvp.yticklabel)
    set(gca,'ytick',1:numel(nvp.yticklabel),'yticklabel',[nvp.yticklabel{:}])
end
a=gca;
a.Interactions = [dataTipInteraction panInteraction];
end

