function [formattedIndices, indx] = formatIndices(indices, matrixSize, h)
    % FORMATINDICES Formats the input indices to eliminate any ':' and turn it into a standard input.
    %
    %   This function formats the input indices to handle special cases like ':' and ':' with a list.
    %
    %   Inputs:
    %       indices    - Cell array of indices.
    %       matrixSize - Size of the matrix (number of rows or columns)
    %       h          - Half the size of the block.
    %
    %   Output:
    %       formattedIndices - Formatted indices.

    blockSize = 2 * h + 1;
    n = matrixSize / blockSize;

    if isequal(indices, {':'})
        formattedIndices = num2cell(1:n);
        formattedIndices = [formattedIndices; repmat({-h:h}, size(formattedIndices))];
        formattedIndices = formattedIndices(:)';
    elseif length(indices) == 2 && isequal(indices{1}, ':')
        pList = indices{2};
        formattedIndices = num2cell(1:n);
        formattedIndices = [formattedIndices; repmat({pList}, size(formattedIndices))];
        formattedIndices = formattedIndices(:)';
    else
        formattedIndices = {};
        for ii = 1:2:length(indices)
            if length(indices{ii}) > 1
                blocks = indices{ii};
                pList = indices{ii+1};
                for block = blocks
                    formattedIndices = [formattedIndices, {block}, {pList}];
                end
            else
                formattedIndices = [formattedIndices, indices(ii:ii+1)];
            end
        end

        for ii = 1:2:length(formattedIndices)
            if isequal(formattedIndices{ii+1}, ':')
                formattedIndices{ii+1} = -h:h;
            end
        end
    end
    if nargout>1
        indx = toIndex(formattedIndices,h);
        indx = unique(indx(:));
    end
    
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