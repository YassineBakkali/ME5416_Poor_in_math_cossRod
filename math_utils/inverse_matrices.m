function out_matrices = inverse_matrices(matrices)
    [rows, cols, n_matrices] = size(matrices);

    % Initialize the output matrix collection with the same size
    out_matrices = zeros(rows, cols, n_matrices);

    % Loop through each matrix along the third dimension
    for n = 1:n_matrices
        % Calculate the inverse of each matrix and store it
        out_matrices(:, :, n) = inv(matrices(:, :, n));
    end