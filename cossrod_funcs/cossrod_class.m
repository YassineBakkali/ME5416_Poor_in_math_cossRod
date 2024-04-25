classdef cossrod_class
    % MyClass A simple example class in MATLAB
    properties
        nElems
        tangents 
        radius 
        lengths 
        volumes
        masses
        rho_medium
        rho_material
        rest_lengths
        dilatations
        dilatation_rate
        rest_voronoi_lengths
        voronoi_dilatations
        internal_stress
        internal_forces
        internal_couple
        internal_torques
        external_forces
        external_torques
        dir_vects
        pos_vects
        vel_vects
        acc_vects
        ang_acc_vects
        sigma
        sigma_rest
        shearMat
        bendMat
        kappa
        rest_kappa
        omega_vects
        mass_second_moment_of_inertia
        inv_mass_second_moment_of_inertia
    end

    methods
        % Constructor method
        function obj = cossrod_class(params)
            if nargin > 0 % Check if any input arguments
                obj.nElems = params.nElements;
                obj.tangents = params.tangents;
                obj.radius = params.radius;
                obj.lengths = params.lengths';
                obj.dilatations = params.dilatation;
                obj.dilatation_rate = params.dilatationRate;
                obj.rest_lengths = params.restLengths;
                obj.voronoi_dilatations = params.voronoiDilatation;
                obj.sigma = params.sigma;
                obj.kappa = params.kappa;
                obj.volumes = params.volume;
                obj.internal_stress = params.internalStress;
                obj.shearMat = params.shearMatrix;
                obj.sigma_rest = params.restSigma;
                obj.internal_forces = params.internalForces;
                obj.internal_couple = params.internalCouple;
                obj.rest_kappa = params.restKappa;
                obj.bendMat = params.bendMatrix;
                obj.pos_vects = params.position;
                obj.vel_vects = params.velocities;
                obj.dir_vects = params.directors;
                obj.omega_vects = params.omegas;
                obj.mass_second_moment_of_inertia = params.massSecondMomentOfInertia;
                obj.inv_mass_second_moment_of_inertia = params.invMassSecondMomentOfInertia;
                obj.internal_torques = params.internalTorques;
                obj.acc_vects = params.accelerations;
                obj.ang_acc_vects = params.angularAccelerations;
                obj.external_forces = params.externalForces;
                obj.masses = params.mass;
                obj.external_torques = params.externalTorques;
                obj.rest_voronoi_lengths = params.restVoronoiLengths;
                obj.rho_material = params.density;              %new
                obj.rho_medium = params.density_medium;         %new

                obj = obj.update_shearStretch(obj.dir_vects, obj.pos_vects,...
                obj.volumes, obj.rest_lengths, obj.rest_voronoi_lengths);

                obj = obj.compute_bending_twist_strains(obj.dir_vects,...
                    obj.rest_voronoi_lengths);
            end
        end

        function obj = compute_geom(obj, pos_vect, volumes)
            %COMPUTE_GEOM Summary of this function goes here :
            %   First, we must compute the difference in position between
            %   all nodes. From there we can get the lengths between nodes
            %   Then use these to get the tangents and updated radi at each node.
            
            pos_diffs =  diff(pos_vect,1,2);
            obj.lengths = vecnorm(pos_diffs) + 1e-14;
            
            obj.tangents = pos_diffs ./ obj.lengths;
            obj.radius = sqrt(volumes ./ obj.lengths ./ pi);
            
        end
            
        
        function obj = compute_dilatations(obj, pos_vect, volumes, rest_lengths, rest_voronoi_lengths)
            %COMPUTE_DILATATIONS Summary of this function goes here
            %   Update the dilatations, voronoi lengths and voronoi dilatations
            
            obj = obj.compute_geom(pos_vect, volumes);
            obj.dilatations = obj.lengths ./ rest_lengths;
            
            voronoi_lengths = (obj.lengths(1:end-1) + obj.lengths(2:end))/ 2;
            obj.voronoi_dilatations = voronoi_lengths ./ rest_voronoi_lengths;

        end

        function obj = update_shearStretch(obj, dir_vects, pos_vect,...
            volumes, rest_lengths, rest_voronoi_lengths)
            %UPDATE_SHEARSTRETCH Summary of this function goes here
            %   Detailed explanation goes here
            obj = obj.compute_dilatations(pos_vect, volumes, rest_lengths, rest_voronoi_lengths);
            
            z_euler = [0; 0; 1];
            obj.sigma = obj.dilatations .* batchMatVec(dir_vects, obj.tangents) - z_euler;
        end

        function obj = update_internal_stresses(obj, dir_vects, pos_vects, volumes,...
                rest_lengths, rest_voronoi_lengths, shearMat, sigma_rest)
            obj = obj.update_shearStretch(dir_vects, pos_vects,...
                volumes, rest_lengths, rest_voronoi_lengths);
            obj.internal_stress = batchMatVec(shearMat, obj.sigma - sigma_rest);
        end

        function obj = compute_internal_forces(obj, dir_vects, pos_vects, volumes,...
                rest_lengths, rest_voronoi_lengths, shearMat, sigma_rest)

            obj = obj.update_internal_stresses(dir_vects, pos_vects, volumes,...
                rest_lengths, rest_voronoi_lengths, shearMat, sigma_rest);

            blocksize = size(obj.internal_stress, 2);
            cosserat_internal_stress = zeros(3, blocksize);
            
        
            % Compute Q^T * n_L / e (transpose directors and multiply by stresses)
            for i = 1:3
                for j = 1:3
                    for k = 1:blocksize
                        cosserat_internal_stress(i, k) = cosserat_internal_stress(i, k) + ...
                            obj.dir_vects(j, i, k) * obj.internal_stress(j, k);
                    end
                end
            end
            
            % Adjust for dilatation, which should be appropriately defined as either scalar or vector
            cosserat_internal_stress = cosserat_internal_stress ./ (obj.dilatations + 1e-14);
            
            % Compute the internal forces using a kernel function
            % Assuming difference_kernel_for_block_structure and ghost_elems_idx are defined appropriately
            obj.internal_forces = two_point_rule(cosserat_internal_stress, obj.nElems);
        end
        

        function obj = compute_bending_twist_strains(obj, director_collection, rest_voronoi_lengths)
            dir_vec = rodrigues_invRot(director_collection);
            obj.kappa = dir_vec ./ rest_voronoi_lengths;
        end

        function obj = compute_internal_bending_twist_strains_from_model(obj)
            obj = obj.compute_bending_twist_strains(obj.dir_vects, obj.rest_voronoi_lengths);
            obj.internal_couple = batchMatVec(obj.bendMat, obj.kappa - obj.rest_kappa);
        end

        function obj = compute_dilatation_rate(obj, pos_vects, vel_vects)
            r_dot_v = sum(pos_vects .* vel_vects, 1);
            r_plus_one_dot_v = sum(pos_vects(:, 2:end) .* vel_vects(:, 1:end-1), 1);
            r_dot_v_plus_one = sum(pos_vects(:, 1:end-1) .* vel_vects(:, 2:end), 1);
            obj.dilatation_rate = (r_dot_v(:,1:end-1) + r_dot_v(:,2:end) - r_dot_v_plus_one - r_plus_one_dot_v) ./ obj.lengths ./ obj.rest_lengths;
        end

        function obj = compute_internal_torques(obj)
            obj = obj.compute_internal_bending_twist_strains_from_model();
            obj = obj.compute_dilatation_rate(obj.pos_vects, obj.vel_vects);
            voronoi_dilatation_inv_cube_cached = 1.0 ./ (obj.voronoi_dilatations .^ 3 + 1e-14);
            bend_twist_couple_2D = two_point_rule(obj.internal_couple .* voronoi_dilatation_inv_cube_cached, []);
            temp = cross(obj.kappa, obj.internal_couple,1);

            bend_twist_couple_3D = trapezoidal_rule(temp .* obj.rest_voronoi_lengths .* voronoi_dilatation_inv_cube_cached);
            shear_stretch_couple = cross(batchMatVec(obj.dir_vects, obj.tangents), obj.internal_stress,1) .* obj.rest_lengths;

            JOmegaUponE = batchMatVec(obj.mass_second_moment_of_inertia, obj.omega_vects) ./ (obj.dilatations + 1e-14);
            lagrangianTransport = cross(JOmegaUponE, obj.omega_vects,1);
        
            % (J \omega_L / e^2) . (de/dt)
            unsteadyDilatation = JOmegaUponE .* obj.dilatation_rate ./ (obj.dilatations + 1e-14);

            obj.internal_torques = bend_twist_couple_2D + bend_twist_couple_3D ...
                + shear_stretch_couple + lagrangianTransport + unsteadyDilatation;
        end

        function obj = compute_internal_forces_and_torques(obj)
            obj = obj.compute_internal_forces(obj.dir_vects, obj.pos_vects, obj.volumes,...
                obj.rest_lengths, obj.rest_voronoi_lengths, obj.shearMat, obj.sigma_rest);
            
            obj = obj.compute_internal_torques();
        end

        function obj = update_acceleration(obj)
            blocksize_acc = size(obj.internal_forces, 2);
            blocksize_alpha = size(obj.internal_torques, 2);
            
            % Calculate accelerations
            for i = 1:3
                for k = 1:blocksize_acc
                    obj.acc_vects(i, k) = ...
                        (obj.internal_forces(i, k) + obj.external_forces(i, k)) / obj.masses(k);
                end
            end
            
            % Initialize alpha_collection to zero
            obj.ang_acc_vects = zeros(size(obj.ang_acc_vects));
            
            % Calculate angular accelerations (alpha)
            for i = 1:3
                for j = 1:3
                    for k = 1:blocksize_alpha
                        obj.ang_acc_vects(i, k) = obj.ang_acc_vects(i, k) + ...
                            (obj.inv_mass_second_moment_of_inertia(i, j, k) * ...
                            (obj.internal_torques(j, k) + obj.external_torques(j, k))) * obj.dilatations(k);
                    end
                end
            end
        end

        function trans_energy = get_translational_energy(obj)
            vel_squared = sum(obj.vel_vects.^2, 1);
            trans_energy = 0.5 * sum(obj.masses .* vel_squared); 
        end

        function rot_energy = get_rotational_energy(obj)
             J_omega_upon_e = batchMatVec(obj.mass_second_moment_of_inertia, obj.omega_vects)...
                 ./ obj.dilatations;
             rot_energy = 0.5 * sum(sum(obj.omega_vects .* J_omega_upon_e, 1));
        end

        function bending_energy = get_bending_energy(obj)
            kappa_diff = obj.kappa - obj.rest_kappa;
            bending_intern_torques = batchMatVec(obj.bendMat, kappa_diff);
            bending_energy = 0.5 * sum(sum(kappa_diff.* bending_intern_torques, 1)...
                .* obj.rest_voronoi_lengths,2);
        end

        function bending_energy = get_shear_energy(obj)
            sigma_diff = obj.sigma - obj.sigma_rest;
            shear_intern_forces = batchMatVec(obj.shearMat, sigma_diff);
            bending_energy = 0.5 * sum(sum(sigma_diff.* shear_intern_forces, 1)...
                .* obj.rest_lengths,2);
        end

    end    
end
