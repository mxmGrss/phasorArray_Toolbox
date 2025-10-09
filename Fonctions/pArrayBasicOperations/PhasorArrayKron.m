function [M] = PhasorArrayKron(Aph,Bph)

arguments
    Aph=[]
    Bph=[]
end

if isa(Aph,'PhasorArray')
    Aph=Aph.Value;
end

if isa(Bph,'PhasorArray')
    Bph=Bph.Value;
end


na=size(Aph,1);
nb=size(Aph,2);

M=[];
for ii=1:na
    Md=[];
    for jj =1:nb
        d=PhasorArrayTimes(Aph(ii,jj,:),Bph);
        Md=[Md d];
    end
    M=[M;Md];
end


end

