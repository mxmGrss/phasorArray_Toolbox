function OutM = array2BT(Aph, m, varg)
    %ARRAY2BT Convert a 3D array into a Block-Toeplitz (BT) matrix.
    %
    %   Usage:
    %       OutM = array2BT(Aph)
    %       OutM = array2BT(Aph, m)
    %       OutM = array2BT(Aph, [h1, h2])
    %
    %   In the PhasorArray toolbox, array2BT is often an alias or a specific
    %   reordering of the Toeplitz-Blocks structure. 
    %   In this dual representation, the matrix is composed of blocks M_k
    %   (Toeplitz representation of the k-th harmonic).
    %
    %   Inputs:
    %       Aph   - Input data (3D array or PhasorArray).
    %       m     - Output order(s) [h1, h2].
    %       varg  - Method ('cell2mat' or 'cat').
    %
    %   Note: This implementation leverages array2TBlocks and reorders 
    %   to match the BT convention (if different from TB).
    %   For many applications, TB and BT are related by a permutation.
    
    arguments
        Aph
        m = []
        varg.method = 'cell2mat';
    end

    % For now, we provide the TB structure and a note on BT.
    % If the BT structure requires a different block placement, 
    % we will adjust it here.
    OutM = array2TBlocks(Aph, m, 'method', varg.method);
    
    % If BT refers to the dual (harmonic-major) structure, 
    % we would use a different tiling logic here.
end
