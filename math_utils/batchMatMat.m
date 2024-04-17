function output_matrix = batchMatMat(first_matrix,second_matrix)
% This function performs batch matrix multiplication for batches of 3x3 matrices.
% Parameters:
% first_matrix_collection : an 3x3xN array of matrices
% second_matrix_collection : an 3x3xN array of matrices
%
% Returns:
% output_matrix : an 3x3xN array of matrix products
%
% Example of usage:
% A = rand(3,3,10); % First set of 3x3 matrices
% B = rand(3,3,10); % Second set of 3x3 matrices
% C = batchMatMul(A, B);
    
    % Determine the size of the third dimension
    blocksize = size(first_matrix, 3);
    
    % Initialize the output matrix collection
    output_matrix = zeros(3, 3, blocksize);
    
    % Loop through each matrix in the batch
    for k = 1:blocksize
        % Matrix multiplication for each pair
        output_matrix(:,:,k) = first_matrix(:,:,k) * second_matrix(:,:,k);
    end
end
