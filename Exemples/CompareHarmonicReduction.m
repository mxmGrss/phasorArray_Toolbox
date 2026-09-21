% Compare coefficient plots and an energy-budget reduction.
% No optional toolbox required. Energy is squared L2, not amplitude error.
amplitudes=[1 .5 .02 .3 .008 .15 .003 .08 .001];
A=PhasorArray(reshape([fliplr(amplitudes(2:end)),amplitudes],1,1,[]));
[B,info]=neglect(A,1e-3,reduceMethod="energy",mode="matrixwise");
fig=figure('Name','Harmonic reduction');
tiles=tiledlayout(fig,2,1,'TileSpacing','compact');
ax=nexttile(tiles);
v=pvalue(A); w=pvalue(B);
stem(ax,0:A.h,abs(reshape(v(:,:,A.h+1:end),1,[])), ...
    'DisplayName','Original'); hold(ax,'on');
stem(ax,0:B.h,abs(reshape(w(:,:,B.h+1:end),1,[])), ...
    'DisplayName','Energy reduction');
title(ax,'Superposed stems'); legend(ax,'show'); grid(ax,'on');
ylabel(ax,'Coefficient magnitude');
bx=nexttile(tiles);
bar(A,B,parent=bx,scale="linear",layout="grouped", ...
    labels=["Original","Energy reduction"]);
title(bx,sprintf('Grouped bars: discarded energy %.4g',info.discardedEnergyFraction));
linkaxes([ax,bx],'xy'); xlim(bx,[-.6,A.h+.6]);
xticks(ax,0:A.h); xlabel(ax,'Harmonic order');
