function out = trapezoidal_rule(vecs)
%TRAPEZOIDAL_RULE Summary of this function goes here
%   Detailed explanation goes here
    [dim, blocksize] = size(vecs);
    tempCollection = zeros(dim, blocksize + 1);
    
    % Apply the trapezoidal rule
    % First and last elements halved
    tempCollection(:, 1) = 0.5 * vecs(:, 1);
    tempCollection(:, end) = 0.5 * vecs(:, end);
    
    % Trapezoidal sums for the internal blocks
    for i = 1:dim
        tempCollection(i, 2:blocksize) = 0.5 * (vecs(i, 1:blocksize-1) + vecs(i, 2:blocksize));
    end
    
    out = tempCollection;
end