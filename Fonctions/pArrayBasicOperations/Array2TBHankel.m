function [HpJ,JHm,Hp,Hm] = Array2TBHankel(Aph,m,varg)
% ARRAY2TBHANKEL Construct Toeplitz Block Hankel matrices from a 3D array.
%
%   ARRAY2TBHANKEL(Aph, m, varg) constructs the Toeplitz Block Hankel (TBH) and J-Hankel
%   matrices from a 3D array containing Fourier coefficients of a periodic matrix function.
%
%   Inputs:
%     Aph     - (3D array) A Fourier coefficient array, typically of size `N×N×(2*nh+1)`.
%     m       - (integer, optional) The truncation order for harmonics.
%                   - Default: The full length of `Aph`.
%     varg    - (optional) Additional parameters:
%                   - 'method' (char): Specifies the assembly method:
%                       - 'cell2mat' (default): Efficient for numeric arrays.
%                       - 'cat': Compatible with symbolic or SDP variables.
%
%   Outputs:
%     HpJ  - ((m+1)N × (m+1)N matrix) Positive J-Hankel matrix.
%     JHm  - ((m+1)N × (m+1)N matrix) Negative J-Hankel matrix.
%     Hp   - ((m+1)N × (m+1)N matrix) Positive Toeplitz Block Hankel matrix.
%     Hm   - ((m+1)N × (m+1)N matrix) Negative Toeplitz Block Hankel matrix.
%
%   Behavior:
%     - The Toeplitz Block Hankel matrices approximate the spectral properties of `Aph`.
%     - If `m` is not provided, the function defaults to the full length of `Aph`.
%
%   Example Usage:
%     % Generate Toeplitz Block Hankel matrices for a given Fourier array
%     Aph = rand(4, 4, 11);  % Fourier coefficients for a 4×4 matrix with 11 harmonics
%     [HpJ, JHm, Hp, Hm] = Array2TBHankel(Aph, 5, 'method', 'cell2mat');
%
%   See also: spArray2TBHankel, TBHankel, BTHankel.
arguments
    Aph
    m=[]
    varg.method {mustBeMember(varg.method,{'cell2mat','cat'})} = 'cell2mat'
end
if isa(Aph,'PhasorArray')
    Aph=Aph.Value;
end



if isa(Aph,"ndsdpvar") || isa(Aph,"sdpvar") || isa(Aph,"sym")
    varg.method='cat';
end
if nargin == 1
    [n1,n2,nhlenbis]=size(Aph);
    nh=(nhlenbis-1)/2;
    m=nh*2;
end

    n1       = size(Aph,1);
    n2       = size(Aph,2);
    nhlenbis = size(Aph,3);
    nh=(nhlenbis-1)/2;
    if nh>2*m+1
        Aph=Aph(:,:,(nh+1+(-2*m-1:2*m+1)));
        nh=2*m+1;
    elseif nh<2*m+1
        switch varg.method
            case 'cell2mat'
        %         Aph=padarray(Aph,[0 0 (m-nh)]);
        %             nh=m
            otherwise
        end
    end



switch varg.method
    case 'cell2mat'
% OutM=zeros(n1*(m+1),n2*(m+1));
cp=cell(n1,n2);
cm=cell(n1,n2);
        for xi = 1:n1
            for yi=1:n2
                ui=[zeros(2*m+1-nh,1) ; squeeze(Aph(xi,yi,:)) ; zeros(2*m+1-nh,1)];
                TAijp = toeplitz(ui((3*m+3):end), ui((3*m+3):-1:(2*m+3)) );
%                 TAij
                cp{xi,yi}=TAijp;
                TAijm = toeplitz(ui((1*m+1):2*m+1), ui((1*m+1):-1:(1)) );
%                 TAij
                cm{xi,yi}=TAijm;
%                 OutM((xi-1)*(m+1)+(1:(m+1)),(yi-1)*(m+1)+(1:(m+1)))=TAij;
            end
        end
        % c
        HpJ=cell2mat(cp);
        JHm=cell2mat(cm);

    case 'cat'
        HpJ=[];
        JHm=[];
        for xi = 1:n1
            HpJ1=[];
            JHm1=[];
            for yi=1:n2
                ui=[zeros(2*m+1-nh,1) ; squeeze(Aph(xi,yi,:)) ; zeros(2*m+1-nh,1)];
                TAijp = toeplitz(ui((3*m+3):end), ui((3*m+3):-1:(2*m+3)) );
                TAijm = toeplitz(ui((1*m+1):2*m+1), ui((1*m+1):-1:(1)) );
                HpJ1=[HpJ1, TAijp];
                JHm1=[JHm1, TAijm];
            end
            HpJ=[HpJ; HpJ1];
            JHm=[JHm; JHm1];
        end
end

J=flip(eye(m+1));
Jp=kron(eye(n2),J);
Jm=kron(eye(n1),J);
Hp=HpJ*Jp;
Hm=Jm*JHm;


end