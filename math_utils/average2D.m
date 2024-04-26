function output_vector = average2D(vector_collection)
    % Compute the size for output_vector
    blocksize = size(vector_collection, 2) - 1;

    % Preallocate output_vector with zeros
    output_vector = zeros(3, blocksize);

    % Calculate averages using vectorized operations
    output_vector = (vector_collection(:, 1:end-1) + vector_collection(:, 2:end)) / 2;
end