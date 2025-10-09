function [pA,Acell,outcat] = BT2array(A,nbloc)
    %BT2ARRAY Convert a Block Toeplitz matrix to a 3D array by extracting diagonal elements
    %   This function takes a Block Toeplitz matrix and extracts the anti-diagonal elements
    %   to convert it into a 3D array. The Block Toeplitz matrix is assumed to have a size defined
    %   by the given block size `nbloc`.
    %
    %   The output consists of three components:
    %   1. A 3D array `pA` containing the phasor values.
    %   2. A cell array `Acell` with the phasor values in cell format.
    %   3. A concatenated version `outcat` where the phasors are combined into a matrix.
    %
    %   Formula:
    %       The input matrix `A` is assumed to be of the form:
    %       [ A0  A-1  A-2 ...  A-(nbloc-1)]
    %       [ A1  A0  A-1  ...  ]
    %       [ A2  A1  A0  ...  ]
    %       ...
    %       [ A(nbloc-2) A(nbloc-3) ... ]
    %       [ Anbloc-1 ...  A0 ]
    %       The function extracts the anti-diagonal elements as a 3D array.
    %
    %   Usage:
    %   [pA, Acell, outcat] = BT2array(A, nbloc)
    %
    %   Inputs:
    %       A      - The Block Toeplitz matrix to be converted.
    %       nbloc  - The block size (dimension of the Toeplitz blocks).
    %
    %   Outputs:
    %       pA     - A 3D array containing the extracted phasors.
    %       Acell  - A cell array containing the phasors.
    %       outcat - A concatenated version of the phasors in matrix form.
    %
    %   See also: num2cell, cell2mat

[nia,nja] = size(A);
nia = nia/nbloc;
nja = nja/nbloc;
pA=zeros(nia,nja,2*nbloc-1);

i0=nbloc;
j0=1;
compteur=(nbloc-1);
while compteur>(-nbloc)
    pA(:,:,compteur+nbloc) = A((i0*nia-nia+1):(i0*nia),(j0*nja-nja+1):(j0*nja));

    if mod(compteur,2) == 0 
        i0=i0-1;
    else
        j0=j0+1;
    end
    compteur=compteur-1;
end
pA;
Acell1= num2cell(pA,[1,2]);
Acell= reshape(Acell1, 1,[]);
outcat=cell2mat(Acell);


end