function [Cph,Bph] = PhasorArrayAdd(Aph)
%PHASORARRAYADD perform padding (and centering) of array of phasor and then
%add their slices to forme A(t) + B(t) in the harmonic domaine. 
%Get 2 3D array and outpout a 3D array
arguments (Repeating)
    Aph
end
n_arg=numel(Aph);

I=cellfun(@(x) isa(x,'PhasorArray'),Aph);
Aph(I)=cellfun(@(x) x.Value,Aph(I),'UniformOutput',false);
I2=cellfun(@(x) isa(x,'sdpvar')||isa(x,'ndsdpvar'),Aph);
I3=cellfun(@(x) isa(x,'sym'),Aph);
% I4=cellfun(@(x) isscalar(x),Aph);


size_arg1=ones(n_arg,1)*1;
size_arg1=num2cell(size_arg1,2)';
size_arg2=ones(n_arg,1)*2;
size_arg2=num2cell(size_arg2,2)';
size_arg3=ones(n_arg,1)*3;
size_arg3=num2cell(size_arg3,2)';

sb1=cellfun(@size,Aph,size_arg1,UniformOutput=false)';
sb2=cellfun(@size,Aph,size_arg2,UniformOutput=false)';
sb3=cellfun(@size,Aph,size_arg3,UniformOutput=false)';

sb=cat(2,sb1,sb2,sb3);

SB=cell2mat(sb);
if range(SB(:,1))~=0 || range(SB(:,2))~=0 
    error("Non compatible dimension")
end



maxh=max(SB(:,3));
toto_1=ones(size(SB,1),1)*[SB(1,1) SB(1,2) maxh];
pad_var=(toto_1-SB)/2; %padding length
padarg=num2cell(pad_var,2)';

if nnz(pad_var)>0 % il y a du padding à faire (sinon ce serai 0)
    if nnz(I2)+nnz(I3)>0 % il a des sdpvar ou des sym dans le tas
        Bph=cell(n_arg,1);
        for kk=1:n_arg
            if SB(kk,3)<maxh % n'a pas assez de phaseurs
                if I2(kk) %Aph{kk} est un sdpvar
                    dummy_ph=sdpvar(SB(1,1),SB(1,2),maxh,'full');
                    dummy_nh0=SB(kk,3); %taille de Aph(kk);
%                     toto(:,:,titi(kk,3)+1:hoho+titi(kk,3))
                    dummy_ph(:,:,pad_var(kk,3)+1:dummy_nh0+pad_var(kk,3))=Aph{kk};
                    dummy_ph(:,:,1:pad_var(kk,3))=0;
                    dummy_ph(:,:,end-pad_var(kk,3)+1:end)=0;
                    Bph{kk}=dummy_ph;
                elseif I3(kk) %Aph{kk} est un sym
                    dummy_ph=sym('toto',[SB(1,1),SB(1,2),maxh]);
                    dummy_nh0=SB(kk,3); %taille de Aph(kk);
%                     toto(:,:,titi(kk,3)+1:hoho+titi(kk,3))
                    dummy_ph(:,:,pad_var(kk,3)+1:dummy_nh0+pad_var(kk,3))=Aph{kk};
                    dummy_ph(:,:,1:pad_var(kk,3))=0;
                    dummy_ph(:,:,end-pad_var(kk,3)+1:end)=0;
                    Bph{kk}=dummy_ph;
                else
                    uu=padarray(Aph{kk},pad_var(kk,:));
                    Bph{kk}=uu;
                end
            else
                Bph{kk}=Aph{kk};
            end
        end
    else
        Bph=cellfun(@padarray,Aph,padarg,UniformOutput=false);
    end
else
    Bph=Aph;
end
    Cph=sum(cat(4,Bph{:}),4);


% if sa(3)>sb(3)
%     Bph=padarray(Bph,[0 0 sa(3)-sb(3)]/2);
% else 
%     Aph=padarray(Aph,[0 0 sb(3)-sa(3)]/2);
% end
% 
% Cph=Aph+Bph;
end

