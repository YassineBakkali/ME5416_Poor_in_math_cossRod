function out = two_point_rule(vecs, ghost_idx)
% This function performs numerical differentiation by computing the difference
% between consecutive elements in a 2D array.
%
% Parameters:
%   arrayCollection - 2D (dim, blocksize) array containing data with 'float' type.
%
% Returns:
%   result - 2D (dim, blocksize-1) array containing data with 'float' type.

    [dim, blocksize] = size(vecs);
    % if ~isempty(ghost_idx)
    %     vecs(:,ghost_idx) = zeros(dim,1);
    % end
    % Get the block size from the second dimension
    blocksize = size(vecs, 2);
    temp_collection = zeros(3, blocksize + 1);  % Using zeros to initialize, equivalent to empty in context
    
    % Setting initial boundary conditions
    temp_collection(1, 1) = vecs(1, 1);
    temp_collection(2, 1) = vecs(2, 1);
    temp_collection(3, 1) = vecs(3, 1);
    
    % Setting final boundary conditions
    temp_collection(1, blocksize + 1) = -vecs(1, blocksize);
    temp_collection(2, blocksize + 1) = -vecs(2, blocksize);
    temp_collection(3, blocksize + 1) = -vecs(3, blocksize);
    
    % Computing differences between successive elements
    for i = 1:3
        for k = 2:blocksize
            temp_collection(i, k) = vecs(i, k) - vecs(i, k - 1);
        end
    end
    out = temp_collection;
end