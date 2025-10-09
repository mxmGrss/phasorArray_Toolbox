function [o2] = PhasorArrayPad(o1,delta_h)
%PhasorArrayPad Pad a PhasorArray with zeros
% PhasorArrayPad(A,h) pads the PhasorArray A with h zeros on each side of the third dimension
% PhasorArrayPad(A,[n1 n2 h]) pads the PhasorArray A with n1 zeros on the first dimension, n2 zeros on the second dimension and h zeros on each side of the third dimension

[n1]=size(o1,1);
[n2]=size(o1,2);
[pre_h_old]=size(o1,3);
h_old=(pre_h_old-1)/2;

if numel(delta_h)==1
    hp = delta_h;
    n1p = 0;
    n2p = 0;
elseif numel(delta_h)==3
    hp = delta_h(3);
    n1p = delta_h(1);
    n2p = delta_h(2);
else 
    error('Invalid input for delta_h in PhasorArrayPad, delta_h must be a scalar or a 3-element vector') 
end

if isphasor(o1)
    output_phas=1;
    o1=o1.value;
else
    output_phas=0;
end
switch class(o1)
    case "double"
        o2=zeros(n1+2*n1p,n2+2*n2p,(2*(h_old+hp)+1));
    case {"ndsdpvar","sdpvar"}
        o2=ndsdpvar(n1+2*n1p,n2+2*n2p,(2*(h_old+hp)+1),'full');
    case "sym"
        o2=sym(zeros(n1+2*n1p,n2+2*n2p,(2*(h_old+hp)+1)));
end
o2(:,:,(1:hp))=0;
o2(n1p+1:n1+n1p,n2p+1:n2+n2p,(hp+1):(end-hp))=o1;
o2(:,:,(end-hp+1):end)=0;

if output_phas
    o2=PhasorArray(o2,reduce=false);
end

end

