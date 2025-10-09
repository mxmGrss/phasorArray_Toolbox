function [Bph] = PhasorUnif(Aph)
%PHASORUNIF take an arbitrary number of phasors arrays, and pad their 3rd
%dimension so they match their number of harmonics
arguments (Repeating)
    Aph
end

I=cellfun(@(x) isa(x,'PhasorArray'),Aph);
Aph(I)=cellfun(@(x) x.value,Aph(I),'UniformOutput',false);

Isdp=cellfun(@(x) isa(x,'ndsdpvar')||isa(x,'sdpvar')||isa(x,'sym'),Aph);

n_arg=numel(Aph);

    size_arg=ones(n_arg,1)*[3];
    size_arg=num2cell(size_arg,2)';
    sb=cellfun(@size,Aph,size_arg,UniformOutput=false)';
    SB(:,3)=cell2mat(sb);
    maxh1=max(SB(:,3));
    toto=ones(size(SB,1),1)*[SB(1,1) SB(1,2) maxh1];
    titi=(toto-SB)/2;
    padarg=num2cell(titi,2)';


if nnz(Isdp)==0 %cas où il n'y a pas de sdpvar dans le tas, c'est facile
    Bph=cellfun(@padarray,Aph,padarg,UniformOutput=false);
else
    for iteri = 1:numel(Aph)
        Aphi=Aph{iteri};
            if isa(Aphi,'ndsdpvar')||isa(Aphi,'sdpvar')
                Aphid=ndsdpvar(size(Aphi,1),size(Aphi,2),maxh1);
                hi=(SB(iteri,3)-1)/2;
                Aphid(:,:,(end+1)/2+(-hi:hi))=Aphi;
                Aphid(:,:,1:(end+1)/2+(-hi-1))=0;
                Aphid(:,:,(end+1)/2+(+hi+1):end)=0;
                Cph{iteri}=Aphid;
            elseif isa(Aphi,'sym')
                Aphid=sym('Aphi',[size(Aphi,1),size(Aphi,2),maxh1]);
                hi=(SB(iteri,3)-1)/2;
                Aphid(:,:,(end+1)/2+(-hi:hi))=Aphi;
                Aphid(:,:,1:(end+1)/2+(-hi-1))=0;
                Aphid(:,:,(end+1)/2+(+hi+1):end)=0;
                Cph{iteri}=Aphid;
            end

%         Bph=cellfun(@padarray,Aph,padarg,UniformOutput=false);
    end
    Cph(~Isdp)=cellfun(@padarray,Aph(~Isdp),padarg(~Isdp),UniformOutput=false);
    Bph=Cph;
end
Dph=Bph;
Bph(I)=cellfun(@(x) PhasorArray(x,reduce=false),Bph(I),'UniformOutput',false);
end
