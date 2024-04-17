function output_vector = batchMatVec(matrix_collection, vector_collection)
    % This function performs batch matrix-vector multiplication
    % matrix_collection is a 3x3xN array of matrices
    % vector_collection is a 3xN array of vectors

    % Determine the block size based on the second dimension of vector_collection
    blocksize = size(vector_collection, 2);
    
    % Pre-allocate output_vector for efficiency
    output_vector = zeros(3, blocksize);
    
    % Perform the batch matrix-vector multiplication using nested loops
    for i = 1:3
        for j = 1:3
            for k = 1:blocksize
                output_vector(i, k) = output_vector(i, k) + ...
                                      matrix_collection(i, j, k) * vector_collection(j, k);
            end
        end
    end
end