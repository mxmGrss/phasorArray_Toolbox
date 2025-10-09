function [O,L,hmlist,Lo,Oo] = Oint(T,ny,n_state,harmlist,gainlist)
%[O,L] = Oint(T,ny,n_state,harmlist,gainlist)
%  Produces matrix O(t) and L(t) of adapted size providing the index of the
%state n_state on y (len(y) = ny), the list of harmonics to be rejected, and a gain (must be a scalar (same gain for all harmonics, or same length as harm_list))
% Syntaxe
% O,L = Oint(T,4,2,[0 1 2],[10 5 2]) produce O of size 5x5 and L 5x4 so that
% dz = O(t) z + Ly is an integrator on composant 2 from y of size 4,  with gain 10 on harmonic 0, 5 on harmonic
% 1, 2 on harm 2.
%
% O,L = Oint(T,4,2,[0 1 2],[10 5 2],1,[0 1],[10 1]) O is now of size 8x8,
% with integrator on 2-component and 1-component of y. On the 2, 0 ,1 and
% 2d harm are integrated, on the 1rst, 0 and 1 are.

arguments
    T
    ny
end
arguments (Repeating)
    n_state
    harmlist
    gainlist
end

N = numel(n_state);

%Motif de base d'un oscillateur de période T/k
matphik = @(k, T) [0 -1; 1 0] * 2 * pi / T * k;

%Matrice d'oscillateurs pour nourrir n oscillateurs pour la liste d'entiers nu, à la freq nu/T
matphi = @(nu, T) arrayfun(matphik, nu, T * ones(numel(nu), 1)', 'UniformOutput', false);


%motif de base pour nourir un oscillateurs
Gbar=[0;1];
% matG =@(gam0,gam) [gam0; kron(ones(numel(gam),1).*(gam(1:end))',Gbar)];
matG = @(gam0, gam) [gam0; kron(gam', Gbar)];

Oo={};
Lo={};
lenL=0;
for ni = 1:N
    harmlist_t = sort(harmlist{ni});
    gainlist_t = gainlist{ni};

    if harmlist_t(1) == 0
        Oo1=0;
        harmlist_t=harmlist_t(2:end);
        nharm=numel(harmlist_t);
        gam0 = gainlist_t(1);
        gamh = gainlist_t(2:end);
        n0 = 1;
    else
        Oo1 = [];
        nharm = numel(harmlist_t);
        n0 = 0;
        gam0 = [];
        gamh = gainlist_t;
    end

    Oo2 = matphi(harmlist_t,T);
    Oo2 = blkdiag(Oo2{:});
    Oo3 = blkdiag(Oo1,Oo2);
    Oo3 = repmat({Oo3},1,numel(n_state{ni}));
    Oo = [Oo,Oo3];


    if isscalar(gainlist_t) %si uniquement un gain est fourni : on l'étend pour l'appliquer à tous les harmoniques
        %     gainlist_t=ones(1,nharm+n0)*gainlist_t;
        gamh = ones(1,nharm)*gainlist_t;
    end

    if numel(gainlist_t)==2
        %     gainlist_t=[gainlist_t(1) , ones(1,nharm)*gainlist_t(2)];
        gam = ones(1,nharm)*gainlist_t(2);
    end

    % Lo{ni}=matG(gainlist_t);
    Lo = [Lo,repmat({matG(gam0,gamh)},1,numel(n_state{ni})) ];
    lenL = lenL+numel(Lo{ni});
end

L = zeros(lenL,ny);
li = 1;
nj = 1;
for ni = 1:N
    for n = n_state{ni}
        L(li:(li+numel(Lo{nj})-1),n) = Lo{nj};
        Looo = L(li:(li+numel(Lo{nj})-1),:);
        li = li+numel(Lo{nj});
        Lo{nj} = Looo;
        nj = nj+1;
    end
end

O = blkdiag(Oo{:});

hmlist = [];
for ii = 1:N
    hmlist =  [ hmlist  repmat(sort(harmlist{ii}),1,numel(n_state{ii}))];
end
end