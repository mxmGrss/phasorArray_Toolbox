function [Mcs,hM] = Phasor2CosSin(Mph,varg)
arguments
    Mph
    varg.only_posk=false
    varg.harmodim=[]
    varg.already_cs=false
end

if ismatrix(Mph)
    if isempty(varg.harmodim)
        warning("Attention : 2D object fourni en tant que phasor, interprété comme un phasor 0 d'une matrice constante, specifier la dimension de stockage des harmonique en cas de vecteur")
    elseif varg.harmodim==2
        Mph=permute(Mph,[1 3 2]);
    elseif varg.harmodim==1
        Mph=permute(Mph,[2 3 1]);
    end
end

if varg.already_cs
    Mcs = Mph;
    hM = (size(Mph,3)-1)/2;
    return
end

n4=size(Mph,4);
n5=size(Mph,5);
for iter_ii=1:n5
    for iter_i=1:n4
        if ~varg.only_posk
            Mph_t=Mph(:,:,(end+1)/2:end,iter_i,iter_ii);
        end
        hM=size(Mph_t,3)-1;
        Mcs(:,:,:,iter_i,iter_ii)=cat(3,-2*imag(flip(Mph_t(:,:,2:end),3)), real(Mph_t(:,:,1)), 2*real(Mph_t(:,:,2:end)));
    end
end
end

