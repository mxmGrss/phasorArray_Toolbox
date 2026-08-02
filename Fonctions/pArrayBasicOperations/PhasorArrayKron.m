function [M] = PhasorArrayKron(Aph,Bph)

arguments
    Aph=[]
    Bph=[]
end

Aph = pvalue(Aph);
Bph = pvalue(Bph);


na=size(Aph,1);
nb=size(Aph,2);

M=[];
for ii=1:na
    Md=[];
    for jj =1:nb
        % Unwrapped before concatenation. PhasorArrayTimes honours "output" on
        % the numeric path but not on the sym one, where it wraps the result
        % anyway, and building the raw block matrix below out of PhasorArray
        % objects fails on a dimension mismatch.
        d=pvalue(PhasorArrayTimes(Aph(ii,jj,:),Bph,[],"output","Array"));
        Md=[Md d];
    end
    M=[M;Md];
end


end

