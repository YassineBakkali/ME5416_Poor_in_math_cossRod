classdef ForwardBackwardAlgo
    properties
        rod
        params
        target
        stepsize
        iter
        endflag
        activation_diff_tol
        muscles
        s_activations
        activations
        prev_activations
        activation_diff_tolerance 

        % Costate parameters
        % Material frame
        internal_force
        internal_couple
        
        % Lab frame
        internal_force_lab 
        internal_couple_lab
        internal_force_derivative
        internal_couple_derivative
    end

    methods
        function obj = ForwardBackwardAlgo(rod, params, muscles, targetPoint)
        obj.rod = rod;
        obj.target = targetPoint;
        obj.params = params;
       

        % Costate definitions

        obj.internal_force = zeros(3, rod.nElems);
        obj.internal_couple = zeros(3, rod.nElems-1);
        obj.internal_force_lab = zeros(3, rod.nElems);
        obj.internal_couple_lab = zeros(3, rod.nElems);
        obj.internal_force_derivative = zeros(3, rod.nElems);
        obj.internal_couple_derivative = zeros(3, rod.nElems);

        % ---
        obj.stepsize = params.stepsize;
        obj.activation_diff_tolerance = params.activation_diff_tol;
        obj.iter = 0;
        obj.endflag = false;

        % Activation and muscles definitions 
        obj.activation_diff_tol = obj.params.activation_diff_tol;
        obj.muscles = muscles;
        obj.s_activations = {};
        obj.activations = {};
        obj.prev_activations = {};
        for idx = 1:length(muscles)
            obj.s_activations{end+1} = muscles{idx}.s_activation;
            obj.activations{end+1} = muscles{idx}.activation;
            obj.prev_activations{end+1} = Inf(size(muscles{idx}.activation));
        end
        end

        function obj = runAlgo(obj, maxIter)
            for idx = 1:maxIter
                obj = obj.update();
                % disp(obj.rod.pos_vects(:,end));
                if obj.endflag
                    fprintf('Finishing the algorithm at iteration %d\n', obj.iter);

                    return;
                end
            end
            fprintf('Finishing the algorithm at max iterations %d\n', obj.iter);
        end

        function obj = update(obj)
            obj.prev_activations= obj.activations;

            obj = obj.find_equi_strains();

            obj = obj.update_from_strain(obj.rod.sigma, obj.rod.kappa);

            obj.target = obj.target.update_costs(obj.rod.pos_vects, obj.rod.dir_vects);

            obj = obj.discrete_cost_gradient_condition();
            
            obj = obj.continuous_cost_gradient_condition();

            obj = obj.costate_backward_evolution();
            [obj, targ_activations] = obj.find_target_activations();
            obj = obj.update_activations(targ_activations);

            obj = obj.check_activations_diff();
            
            obj.iter = obj.iter + 1;

        end

        function obj = check_activations_diff(obj)
            
            normVal = 0;  % Initialize norm to zero
            numActivations = numel(obj.activations);  % Determine the number of activations
            
            for i = 1:numActivations
                diff = obj.activations{i} - obj.prev_activations{i};
                normVal = normVal + sum(diff.^2) / numel(obj.activations{i});
            end
            
            % Average the accumulated norm values
            normVal = normVal / numActivations;
            
            % Check against the tolerance and return the result
            obj.endflag = normVal < obj.activation_diff_tolerance;
        end

        function obj = update_activations(obj, target_activations)
            
            numActivations = numel(obj.activations);  % Number of activations
            
            for i = 1:numActivations
                % Update each activation according to the target
                obj.activations{i} = obj.activations{i} - obj.stepsize * (obj.activations{i} - target_activations{i});
                
                % Apply clipping to the updated activations to ensure they are within [0, 1]
                obj.activations{i} = max(min(obj.activations{i}, 1), 0);
            end
        end

        function [obj, target_activations] = find_target_activations(obj)
            numMuscles = numel(obj.muscles); 
            target_activations = cell(1, numMuscles);  
            
            for i = 1:numMuscles
                obj.muscles{i} = obj.muscles{i}.set_activation(ones(size(obj.muscles{i}.activation)));  % Apply uniform activation
                obj.muscles{i} = obj.muscles{i}.update_muscle_force_couple(obj.rod);  % Update muscle state based on static rod
                
                % Calculate the target activation and store it
                [target_activations{i}, obj] = obj.calculate_target_activations(...
                    obj.muscles{i}.internal_forces, ...
                    obj.muscles{i}.internal_couples);
            end
        end

        function [target_activation, obj] = calculate_target_activations(obj, mscl_internal_force, mscl_internal_couple)
            blocksize = size(obj.internal_force, 2);
            target_activation = zeros(1, blocksize);
            
            % Perform matrix-vector multiplication with the inverse matrices
            temp_shear = batchMatVec(inverse_matrices(obj.rod.shearMat./reshape(obj.rod.dilatations,[1,1,obj.rod.nElems])), mscl_internal_force);
            temp_kappa = batchMatVec(inverse_matrices(obj.rod.bendMat), mscl_internal_couple) .* (obj.rod.voronoi_dilatations);
            
            temp_force_innerproduct = zeros(1, blocksize);
            temp_couple_innerproduct = zeros(1, blocksize + 1);
            
            % Compute the inner products for force
            for k = 1:blocksize
                for i = 1:3
                    temp_force_innerproduct(k) = temp_force_innerproduct(k) + ...
                        obj.internal_force(i, k) * temp_shear(i, k);
                end
            end
            
            % Compute the inner products for couple
            for k = 1:blocksize-1
                for i = 1:3
                    temp_couple_innerproduct(k+1) = temp_couple_innerproduct(k+1) + ...
                        obj.internal_couple(i, k) * temp_kappa(i, k);
                end
            end
            
            % Apply boundary conditions for the couple inner product
            temp_couple_innerproduct(1) = 2 * temp_couple_innerproduct(2) - temp_couple_innerproduct(3);
            temp_couple_innerproduct(end) = 2 * temp_couple_innerproduct(end-1) - temp_couple_innerproduct(end-2);
            
            % Calculate the target activation
            for k = 1:blocksize
                target_activation(k) = -temp_force_innerproduct(k) - 0.5 * ...
                    (temp_couple_innerproduct(k) + temp_couple_innerproduct(k+1));
            end
        end

        function obj = discrete_cost_gradient_condition(obj)
            obj.internal_force_lab(:, end) = ...
                obj.target.dis_grad_pose(:, end);

            obj.internal_couple_lab(:, end) = ...
                -obj.target.dis_grad_dir(:,end);
        end
        
        function obj = continuous_cost_gradient_condition(obj)
            obj.internal_force_derivative(:, :) = ...
                obj.target.cont_grad_pose;
            obj.internal_couple_derivative(:, :) = ...
                obj.target.cont_grad_dir;
        end

        function obj = find_equi_strains(obj)
            [obj, temp_muscle_forces, temp_muscle_couples] = obj.calc_muscle_forceCouple();
            obj.rod.kappa = - batchMatVec(inverse_matrices(obj.rod.bendMat), temp_muscle_couples).*(obj.rod.voronoi_dilatations.^3);

            obj.rod.sigma = - batchMatVec(inverse_matrices(obj.rod.shearMat./reshape(obj.rod.dilatations,[1,1,obj.rod.nElems])),...
                temp_muscle_forces);
            
        end

        function obj = costate_backward_evolution(obj)
    
            blocksize = length(obj.rod.rest_lengths);
            internal_force_lab_frame = zeros(3, blocksize);
            internal_couple_lab_frame = zeros(3, blocksize);
            
            % Calculate shear and dilatation
            shear = obj.rod.sigma;
            shear(3, :) = obj.rod.sigma(3, :) + 1;
            temp = shear;
            temp(3, :) = temp(3, :) + 1;

            dilatation = vecnorm(temp);
            
            internal_force_lab_frame(:, end) = obj.internal_force_lab(:, end);
            
            % Average force derivatives
            force_derivative = average2D(obj.internal_force_derivative);
            
            % Calculate internal force lab frame
            for k = 0:blocksize-2
                internal_force_lab_frame(:, end-k-1) = internal_force_lab_frame(:, end-k) - ...
                    (force_derivative(:, end-k) .* 0.5 .* (obj.rod.rest_lengths(end-k) .* dilatation(end-k) + obj.rod.rest_lengths(end-k-1) .* dilatation(end-k-1)));
            end
            obj.internal_force = batchMatVec(obj.rod.dir_vects, internal_force_lab_frame);
            
            % Set initial conditions for internal couple
            internal_couple_lab_frame(:, end) = obj.internal_couple_lab(:, end);
            
            % Calculate internal couple derivatives and update lab frame couples
            obj.internal_couple_derivative = obj.internal_couple_derivative - cross(...
                material2labframe(obj.rod.dir_vects, shear), internal_force_lab_frame,1);
            couple_derivative = average2D(obj.internal_couple_derivative);
            
            for k = 0:blocksize-2
                internal_couple_lab_frame(:, end-k-1) = internal_couple_lab_frame(:, end-k) - ...
                    (couple_derivative(:, end-k) * 0.5 * (obj.rod.rest_lengths(end-k) * dilatation(end-k) + obj.rod.rest_lengths(end-k-1) * dilatation(end-k-1)));
            end
            obj.internal_couple = average2D(batchMatVec(obj.rod.dir_vects, internal_couple_lab_frame));
        end

        function [obj, muscle_forces, muscle_couples] = calc_muscle_forceCouple(obj)
            muscle_forces = zeros(size(obj.rod.sigma));
            muscle_couples = zeros(size(obj.rod.kappa));

            for idx = 1:length(obj.muscles)
                activation = obj.activations{idx};
    
                obj.muscles{idx} = obj.muscles{idx}.set_activation(activation);
                obj.muscles{idx} = obj.muscles{idx}.update_muscle_force_couple(obj.rod);
                muscle_forces = muscle_forces + obj.muscles{idx}.internal_forces; 
                muscle_couples = muscle_couples + obj.muscles{idx}.internal_couples;
            end
        end

        function obj = update_from_strain(obj, sigma, kappa)
            obj.rod.sigma = sigma;
            obj.rod.kappa = kappa;

            obj = obj.static_pose_evolution();
            obj.rod = obj.rod.compute_dilatations(obj.rod.pos_vects, obj.rod.volumes,...
                obj.rod.rest_lengths, obj.rod.rest_voronoi_lengths);
            
        end

        function obj = static_pose_evolution(obj)
            shear = obj.rod.sigma;
            shear(3, :) = obj.rod.sigma(3, :) + 1;
            % disp(obj.rod.dir_vects(:,:,end))
            for k = 1:size(obj.rod.rest_lengths, 2)-1
                obj.rod.pos_vects(:, k:k+1) = obj.next_pos(obj.rod.dir_vects(:,:,k), shear(:, k) .* obj.rod.rest_lengths(k), ...
                    obj.rod.pos_vects(:, k:k+1));
                obj.rod.dir_vects(:, :, k:k+1) = obj.next_dir(obj.rod.kappa(:, k) .* obj.rod.rest_lengths(k), ...
                      obj.rod.dir_vects(:, :, k:k+1));
            end
            % disp(obj.rod.dir_vects(:,:,end))
            obj.rod.pos_vects(:,end-1:end) = obj.next_pos(obj.rod.dir_vects(:,:,end), ...
                        shear(:, end) .* obj.rod.rest_lengths(end), ...
                        obj.rod.pos_vects(:,end-1:end));

        end

        function out_pos = next_pos(obj, dir, delta, pos)
            out_pos = pos;
            out_pos(:, 2) = pos(:, 1);
        
            
            for i = 1:3
                for j = 1:3
                    out_pos(i, 2) = out_pos(i, 2) + dir(j, i) .* delta(j);
                end
            end
        end

        function out_dir = next_dir(obj, rot, dir)
            rot_mat = rodrigues_roth(reshape(rot, [3, 1]), 1); 
            rot_mat = squeeze(rot_mat(:,:,1));
            out_dir = dir;
            temp_dir = zeros(size(dir(:,:,1))); 
        
            for i = 1:3
                for j = 1:3
                    for k = 1:3
                        temp_dir(i, j) = temp_dir(i, j) + rot_mat(i, k) * dir(k, j, 1);
                    end
                end
            end
            out_dir(:,:,2) = temp_dir;
        end
    end
end
