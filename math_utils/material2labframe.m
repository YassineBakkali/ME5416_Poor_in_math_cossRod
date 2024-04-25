function labFrameVects = material2labframe(dir_vects, vects)
    % Get the size of the second dimension of vectors
    blocksize = size(vects, 2);
    
    % Initialize lab_frame_vectors with zeros
    labFrameVects = zeros(3, blocksize);
    
    % Perform matrix multiplication for each slice n
    for n = 1:blocksize
        labFrameVects(:, n) = dir_vects(:, :, n) * vects(:, n);
    end
end