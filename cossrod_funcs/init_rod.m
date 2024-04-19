function params = init_rod(...
    nElements, direction, normal, baseLength, baseRadius, density, density_medium, youngsModulus, ...
    rodOriginPosition, shearModulus, position, directors, restSigma, restKappa)
    
    % Constants and initial sanity checks
    assert(nElements > 1, 'Invalid number of elements.');
    assert(baseLength > 1e-14, 'Base length is too small.');
    assert(norm(normal) >  1e-14, 'Normal vector is too small.');
    assert(norm(direction) >  1e-14, 'Direction vector is too small.');

    % Define number of nodes and Voronoi elements
    nNodes = nElements + 1;
    nVoronoiElements = nElements - 1;

    % Initialize position
    if isempty(position)
        position = zeros(3, nNodes);
        % Straight open rod
        start = rodOriginPosition;
        endPos = start + direction * baseLength;
        for i = 1:3
            position(i, :) = linspace(start(i), endPos(i), nNodes);
        end
    end

    % Calculate rest lengths and tangents
    positionForDifference = position;

    positionDiff = diff(positionForDifference, 1, 2);
    restLengths = vecnorm(positionDiff);
    tangents = positionDiff ./ restLengths;
    normal = normal / norm(normal);

    % Initialize directors
    if isempty(directors)
        % Set the directors matrix with dimensions 3x3xnElements, assuming 3D space
        directors = zeros(3, 3, nElements);
        
        % Construct directors using tangents and normal
        normalCollection = repmat(normal, 1, nElements);  % Replicating the normal vector across columns
        disp(normalCollection)
        % Check if rod normal and rod tangent are perpendicular to each other
        % Compute dot product across all elements to ensure perpendicularity
        dotProducts = sum(normalCollection .* tangents, 1);  % element-wise multiplication and sum over rows
        assert(all(abs(dotProducts) < 1e-14), ...
            'Rod normal and tangent are not perpendicular to each other!');
        
        % Assign normal, cross product of tangents and normal, and tangents to the directors matrix
        directors(1, :, :) = normalCollection;  % Assign normals to first row of each director slice
        for k = 1:nElements
            directors(2, :, k) = cross(tangents(:, k), normalCollection(:, k));  % Cross product
            directors(3, :, k) = tangents(:, k);  % Tangents
        end
    end

    % Set radius and density arrays
    radius = repmat(baseRadius, nElements, 1);
    assert(all(radius > 1e-14), 'Radius must be greater than 0.');
    densityArray = repmat(density, nNodes, 1);
    assert(all(densityArray > 1e-14), 'Density must be greater than 0.');
    densityMediumArray = repmat(density_medium, nNodes, 1);
    assert(all(densityMediumArray > 1e-14), 'Density must be greater than 0.');

    A0 = pi * radius.^2;
    I0_1 = A0.^2 / (4 * pi);
    I0_2 = I0_1;
    I0_3 = 2 * I0_2;
    I0 = [I0_1, I0_2, I0_3];
    
    % Mass second moment of inertia for disk cross-section
    massSecondMomentOfInertia = zeros(3, 3, nElements);
    for i = 1:nElements
        massSecondMomentOfInertia(:, :, i) = diag([I0(i, 1), I0(i, 2), I0(i, 3)] .* density * restLengths(i));
    end
    
    % Inverse of second moment of inertia
    invMassSecondMomentOfInertia = zeros(3, 3, nElements);
    for i = 1:nElements
        invMassSecondMomentOfInertia(:, :, i) = inv(massSecondMomentOfInertia(:, :, i));
    end

    % Shear modulus computation, default handling
    if isempty(shearModulus)
        poissonRatio = 0.5; % Default Poisson's ratio if shear modulus is not provided
        shearModulus = youngsModulus / (2 * (1 + poissonRatio));
    end
    
    % Shear/Stretch matrix computation
    alphaC = 27 / 28; % Based on Timoshenko correction factor
    shearMatrix = zeros(3, 3, nElements);
    for i = 1:nElements
        shearMatrix(:, :, i) = diag([alphaC * shearModulus * A0(i), alphaC * shearModulus * A0(i), youngsModulus * A0(i)]);
    end

    % Bend/Twist matrix computation
    bendMatrix = zeros(3, 3, nVoronoiElements + 1);
    for i = 1:nElements
        bendMatrix(:, :, i) = diag([youngsModulus * I0_1(i), youngsModulus * I0_2(i), shearModulus * I0_3(i)]);
    end

    % Compute volume of elements
    radius = radius';
    volume = pi * radius.^2 .* restLengths;
    
    % Compute mass of elements
    mass = zeros(1, nNodes); % Initialize mass array for all nodes

    mass(1:end-1) = mass(1:end-1) + 0.5 * density .* volume; % Distribute mass to the beginning node of each element
    mass(2:end) = mass(2:end) + 0.5 * density .* volume;    % Distribute mass to the ending node of each element

    
    % Generate rest sigma and rest kappa, use user input if defined
    % Set rest strains and curvature to be zero at start
    % Initialize rest_sigma and rest_kappa based on provided input or default to zeros
    if isempty(restSigma)
        restSigma = zeros(3, nElements); % Assuming 3D (change dimensions as needed)
    else
        % Optionally validate shape if needed
        assert(isequal(size(restSigma), [3, nElements]), 'restSigma must have the size [3, nElements]');
    end
    
    if isempty(restKappa)
        restKappa = zeros(3, nVoronoiElements); % Assuming 3D (change dimensions as needed)
    else
        % Optionally validate shape if needed
        assert(isequal(size(restKappa), [3, nVoronoiElements]), 'restKappa must have the size [3, nVoronoiElements]');
    end

    % Calculate rest lengths for Voronoi
    restVoronoiLengths = 0.5 * (restLengths(1:end-1) + restLengths(2:end));

    % Initialize velocities and angular velocities
    velocities = zeros(3, nNodes);
    omegas = zeros(3, nElements);
    accelerations = zeros(3, nNodes);
    angularAccelerations = zeros(3, nElements);

    % Initialize internal and external forces and torques
    internalForces = zeros(3, nNodes);
    internalTorques = zeros(3, nElements);
    externalForces = zeros(3, nNodes);
    externalTorques = zeros(3, nElements);

    % Initialize other properties
    lengths = zeros(1, nElements);
    tangents = zeros(3, nElements);
    dilatation = zeros(1, nElements);
    voronoiDilatation = zeros(1, nVoronoiElements);
    dilatationRate = zeros(1, nElements);

    sigma = zeros(3, nElements);
    kappa = zeros(3, nVoronoiElements);
    internalStress = zeros(3, nElements);
    internalCouple = zeros(3, nVoronoiElements);

    % Organizing all the results into a struct 'params'
    params = struct(...
        'nElements', nElements, ...
        'position', position, ...
        'velocities', velocities, ...
        'omegas', omegas, ...
        'accelerations', accelerations, ...
        'angularAccelerations', angularAccelerations, ...
        'directors', directors, ...
        'radius', radius, ...
        'massSecondMomentOfInertia', massSecondMomentOfInertia, ...
        'invMassSecondMomentOfInertia', invMassSecondMomentOfInertia, ...
        'shearMatrix', shearMatrix, ...
        'bendMatrix', bendMatrix, ...
        'density', densityArray, ...
        'density_medium', densityMediumArray, ...
        'volume', volume, ...
        'mass', mass, ...
        'internalForces', internalForces, ...
        'internalTorques', internalTorques, ...
        'externalForces', externalForces, ...
        'externalTorques', externalTorques, ...
        'lengths', lengths, ...
        'restLengths', restLengths, ...
        'tangents', tangents, ...
        'dilatation', dilatation, ...
        'dilatationRate', dilatationRate, ...
        'voronoiDilatation', voronoiDilatation, ...
        'restVoronoiLengths', restVoronoiLengths, ...
        'sigma', sigma, ...
        'kappa', kappa, ...
        'restSigma', restSigma, ...
        'restKappa', restKappa, ...
        'internalStress', internalStress, ...
        'internalCouple', internalCouple ...
    );

end