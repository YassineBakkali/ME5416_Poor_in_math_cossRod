function force_weight = force_length_weight(muscle_length, f_l_coefficients)
    if nargin < 2
        % Default polynomial coefficients if not provided
        f_l_coefficients = [-6.44, 18.01, -13.64, 3.06];
    end
    
    % Get the number of coefficients (degree + 1 of the polynomial)
    degree = numel(f_l_coefficients);

    % Initialize the force_weight array with zeros
    blocksize = numel(muscle_length);
    force_weight = zeros(1, blocksize);

    % Calculate the polynomial value for each muscle length
    for i = 1:blocksize
        for power = 0:degree-1
            force_weight(i) = force_weight(i) + f_l_coefficients(power + 1) * (muscle_length(i) ^ power);
        end

        % Apply constraints
        if force_weight(i) < 0 || muscle_length(i) > 2
            force_weight(i) = 0;
        end
    end
end