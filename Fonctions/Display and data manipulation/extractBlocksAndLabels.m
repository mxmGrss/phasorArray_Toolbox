function [extractedMatrix, reducedRowLabels, reducedColLabels, Idx, Idy] = extractBlocksAndLabels(A, rowIndices, colIndices, h, rowLabels, colLabels)
% EXTRACTBLOCKSANDLABELS Extracts specified blocks and their labels from a block matrix A of N x M blocks of size (2h+1)*(2*h+1) each.
%
%   [extractedMatrix, reducedRowLabels, reducedColLabels, Idx, Idy] = extractBlocksAndLabels(A, rowIndices, colIndices, h, rowLabels, colLabels)
%
%   This function extracts specified blocks and index inside blocks from a matrix A based on the provided
%   row and column indices. It also extracts or generates corresponding labels for
%   the rows and columns of the extracted blocks.
%   Each block of A is assumed square of size (2h+1) x (2h+1).
%   A is assumed to be a block matrix of size N x M, where N and M are the number of blocks in the row and column directions.
%   Each block is indexed by a row-block index and a column-block index nx, ny in the range [1, N] and [1, M],
%   and each of their component is sub indexed by p, q in the range [-h, h].
%   Any component of A is then indexed by the quadruplet (nx, p, ny, q).
%
%   This function can be usefull to extract specific component from a Toeplitz-block matrix.
%
%   Inputs:
%       A           - The block matrix of size (nx*(2h+1)) x (ny*(2h+1)).
%       rowIndices  - Cell array in the form {nxi, [list of p], nxii, [list of p for nxii], ...}
%                     where nxi is the block row index and [list of p] are the indices within the block.
%                       the indices p are in the range [-h, h].
%                     Special cases:
%                       - {':'} takes all rows.
%                       - {':', [list of p]} uses the same list for every block.
%                       - {nxi, ':', nxii, [list of p]} takes all rows for nxi, but selective for nxii.
%       colIndices  - Cell array in the form {nyi, [list of q], ...}
%                     where nyi is the block column index and [list of q] are the indices within the block.
%                       the indices q are in the range [-h, h].
%                     Special cases:
%                       - {':'} takes all columns.
%                       - {':', [list of q]} uses the same list for every block.
%                       - {nyi, ':', nyii, [list of q]} takes all columns for nyi, but selective for nyii.
%       h           - Half the size of the block (block size is (2h+1) x (2h+1)).
%       rowLabels   - (Optional) Cell array of size ny*(2h+1) containing row labels. If empty, labels will be generated.
%       colLabels   - (Optional) Cell array of size nx*(2h+1) containing column labels. If empty, labels will be generated.
%
%   Outputs:
%       extractedMatrix  - The extracted matrix containing the specified blocks.
%       reducedRowLabels - Cell array containing the labels for the rows of the extracted matrix.
%       reducedColLabels - Cell array containing the labels for the columns of the extracted matrix.
%       Idx              - Row indices used to extract the blocks from A.
%       Idy              - Column indices used to extract the blocks from A.
%
%   Example:
%       h = 2;
%       A = reshape(1:(3 * 3 * (2 * h + 1)^2), 3 * (2 * h + 1), 3 * (2 * h + 1));
%       rowIndices = {1, ':', 2, [-1, 1]};
%       colIndices = {1, [0, 1], 3, [-1, 0]};
%       [extractedMatrix, reducedRowLabels, reducedColLabels, Idx, Idy] = extractBlocksAndLabels(A, rowIndices, colIndices, h);
%       disp('Extracted Matrix:');
%       disp(extractedMatrix);
%       disp('Reduced Row Labels:');
%       disp(reducedRowLabels);
%       disp('Reduced Col Labels:');
%       disp(reducedColLabels);
%       disp('Row Indices:');
%       disp(Idx);
%       disp('Column Indices:');
%       disp(Idy);

% Validate input arguments
arguments
    A (:,:) double
    rowIndices {mustBeA(rowIndices, 'cell')}
    colIndices {mustBeA(colIndices, 'cell')}
    h (1,1) double {mustBePositive, mustBeInteger}
    rowLabels (:,:) cell = {}
    colLabels (:,:) cell = {}
end

% Calculate the size of each block
blockSize = 2 * h + 1;

% Initialize empty arrays to store the indices
Idx = [];
Idy = [];

% Format row and column indices
rowIndices = formatIndices(rowIndices, size(A, 1), h);
colIndices = formatIndices(colIndices, size(A, 2), h);

% Generate default labels if they are empty
if isempty(rowLabels)
    labels = arrayfun(@(p) arrayfun(@(q) {sprintf('y_{%d,%d}', p, q)},-h:h), 1:size(A, 1)/blockSize , 'UniformOutput', false);
    rowLabels = [labels{:}];
end
if isempty(colLabels)
    labels = arrayfun(@(p) arrayfun(@(q) {sprintf('u_{%d,%d}', p, q)}, -h:h), 1:size(A, 2)/blockSize, 'UniformOutput', false);
    colLabels = [labels{:}];
end

Idx = toIndex(rowIndices, h);
Idy = toIndex(colIndices, h);


% Extract the submatrix from A
extractedMatrix = A(Idx, Idy);

% Select the labels based on the indices
reducedRowLabels = rowLabels(Idx);
reducedColLabels = colLabels(Idy);
end


function Idx = toIndex(indxCell, h)
Idx = [];
% Generate row indices
blockSize = 2 * h + 1;
for i = 1:2:length(indxCell)
    blockRow = indxCell{i};
    pList = indxCell{i+1};

    % Map the p indices from [-h, h] to [1, 2h+1]
    pMapped = pList + h + 1;

    % Calculate the starting indices for the block
    rowStart = (blockRow - 1) * blockSize + 1;

    % Store the row indices
    tempIdx = rowStart + pMapped - 1;
    Idx = [Idx; tempIdx'];
end
end