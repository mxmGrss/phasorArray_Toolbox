function [Cph,Bph] = PhasorArrayBilAdd(Aph,Bph)
%PHASORARRAYADD perform padding (and centering) of array of phasor and then
%add their slices to forme A(t) + B(t) in the harmonic domaine. 
%Get 2 3D array and outpout a 3D array
arguments 
    Aph
    Bph
end
n_arg=numel(Aph);

if isa(Aph,'PhasorArray')
    Aph=Aph.value;
end

if isa(Bph,'PhasorArray')
    Bph=Bph.value;
end

I=cellfun(@(x) isa(x,'PhasorArray'),Aph);
Aph(I)=cellfun(@(x) x.Value,Aph(I),'UniformOutput',false);


size_arg=ones(n_arg,1)*[1 2 3];
size_arg=num2cell(size_arg,2)';

sb=cellfun(@size,Aph,size_arg,UniformOutput=false)';



SB=cell2mat(sb);

if range(SB(:,1))~=0 || range(SB(:,2))~=0 
    error("Non compatible dimension")
end



maxh=max(SB(:,3));
toto=ones(size(SB,1),1)*[SB(1,1) SB(1,2) maxh];
titi=(toto-SB)/2;
padarg=num2cell(titi,2)';
Bph=cellfun(@padarray,Aph,padarg,UniformOutput=false);
Cph=sum(cat(4,Bph{:}),4);


% if sa(3)>sb(3)
%     Bph=padarray(Bph,[0 0 sa(3)-sb(3)]/2);
% else 
%     Aph=padarray(Aph,[0 0 sb(3)-sa(3)]/2);
% end
% 
% Cph=Aph+Bph;
end

