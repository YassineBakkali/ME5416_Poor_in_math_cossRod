function [rot_mat] = rodrigues_roth(axis_vecs, scale)
%RODRIGUES_ROTH Summary of this function goes here
%   Detailed explanation goes here
    % Get the number of blocks (vectors in the collection)
    blocksize = size(axis_vecs, 2);

    % Calculate the magnitude of each vector
    theta = sqrt(sum(axis_vecs.^2, 1));

    % Avoid division by zero
    small_number = 1e-14;
    theta(theta < small_number) = small_number;

    % Normalize the vectors
    norm_vecs = axis_vecs ./ theta;

    % Scale the angles
    theta = theta * scale;

    % Precompute sine and cosine terms
    u_prefix = sin(theta);
    u_sq_prefix = 1 - cos(theta);

    % Compute cross-product matrix components for each vector
    v0 = norm_vecs(1, :);
    v1 = norm_vecs(2, :);
    v2 = norm_vecs(3, :);

    % Allocate rotation matrix array
    rot_mat = zeros(3, 3, blocksize);

    % Fill the diagonal elements
    rot_mat(1, 1, :) = 1 - u_sq_prefix .* (v1.^2 + v2.^2);
    rot_mat(2, 2, :) = 1 - u_sq_prefix .* (v0.^2 + v2.^2);
    rot_mat(3, 3, :) = 1 - u_sq_prefix .* (v0.^2 + v1.^2);

    % Fill the off-diagonal elements
    rot_mat(1, 2, :) = u_prefix .* v2 + u_sq_prefix .* v0 .* v1;
    rot_mat(2, 1, :) = -u_prefix .* v2 + u_sq_prefix .* v0 .* v1;
    rot_mat(1, 3, :) = -u_prefix .* v1 + u_sq_prefix .* v0 .* v2;
    rot_mat(3, 1, :) = u_prefix .* v1 + u_sq_prefix .* v0 .* v2;
    rot_mat(2, 3, :) = u_prefix .* v0 + u_sq_prefix .* v1 .* v2;
    rot_mat(3, 2, :) = -u_prefix .* v0 + u_sq_prefix .* v1 .* v2;

    % blocksize = size(axis_vecs, 2);
    % rot_mat = zeros(3, 3, blocksize);
    % 
    % for k = 1:blocksize
    %     v0 = axis_vecs(1, k);
    %     v1 = axis_vecs(2, k);
    %     v2 = axis_vecs(3, k);
    % 
    %     % Calculate the magnitude of the vector
    %     theta = sqrt(v0^2 + v1^2 + v2^2);
    % 
    %     % Normalize the vector components
    %     v0 = v0 / (theta + 1e-14);
    %     v1 = v1 / (theta + 1e-14);
    %     v2 = v2 / (theta + 1e-14);
    % 
    %     % Scale the angle by the scale factor
    %     theta = theta * scale;
    %     u_prefix = sin(theta);
    %     u_sq_prefix = 1 - cos(theta);
    % 
    %     % Fill the rotation matrix for this vector
    %     rot_mat(1, 1, k) = 1 - u_sq_prefix * (v1^2 + v2^2);
    %     rot_mat(2, 2, k) = 1 - u_sq_prefix * (v0^2 + v2^2);
    %     rot_mat(3, 3, k) = 1 - u_sq_prefix * (v0^2 + v1^2);
    % 
    %     rot_mat(1, 2, k) = u_prefix * v2 + u_sq_prefix * v0 * v1;
    %     rot_mat(2, 1, k) = -u_prefix * v2 + u_sq_prefix * v0 * v1;
    %     rot_mat(1, 3, k) = -u_prefix * v1 + u_sq_prefix * v0 * v2;
    %     rot_mat(3, 1, k) = u_prefix * v1 + u_sq_prefix * v0 * v2;
    %     rot_mat(2, 3, k) = u_prefix * v0 + u_sq_prefix * v1 * v2;
    %     rot_mat(3, 2, k) = -u_prefix * v0 + u_sq_prefix * v1 * v2;
    % end
end

