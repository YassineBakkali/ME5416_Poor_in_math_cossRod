classdef simulator_class 
    %SIMULATOR_CLASS Summary of this class goes here
    %   Detailed explanation goes here

    properties
        rod
        kine_states
        dyna_states
        kinematic_rates
    end

    methods
        function obj = simulator_class(cosserat_rod)
        %SIMULATOR_CLASS Construct an instance of this class
        %   Detailed explanation goes here
        obj.rod = cosserat_rod;
        obj.kine_states = kine_state(obj.rod.pos_vects, obj.rod.dir_vects);
        obj.dyna_states = dyna_state(zeros(2,3,obj.rod.nElems+1), zeros(2,3,obj.rod.nElems+1), obj.rod.vel_vects, obj.rod.omega_vects);
        obj.kinematic_rates = obj.dyna_states.kinematic_rates();
        end

        function obj = update_internal_forces_and_torques(obj)
            obj.rod.compute_internal_forces_and_torques();
        end

        function dynamic_rates = dynamic_rates(obj, prefac)
            obj.rod.update_acceleration();
            obj.dyna_states.dvdt_dwdt_vecs = permute(cat(3,obj.rod.acc_vects(:,1:obj.rod.nElems), obj.rod.ang_acc_vects), [3, 1, 2]);
            dynamic_rates = obj.dyna_states.dynamic_rates(prefac);
        end

        function obj = integrate(obj, final_time, n_steps)
        %METHOD1 Summary of this method goes here
        %   Detailed explanation goes here
        assert(final_time > 0, 'Final time cannot be negative');
        assert(n_steps > 0, 'Number of integration steps is negative!');
        
        dt = final_time / n_steps;
        time = 0.0;

        for i = 1:n_steps
            time = obj.sim_step(time, dt);
        end

        fprintf(' Simulation time is %.6f\n', time);
        
        end

        function dt = prefac(obj, dt)
            dt = 0.5 * dt;
        end

        function obj = kinematic_step(obj, time, dt)
            prefac = obj.prefac(dt);
            [out_pos, out_dir] = update_kinematics(prefac, obj.kine_states.pos_vecs, ...
                obj.kine_states.dir_vecs, obj.rod.vel_vects, obj.rod.omega_vects);
            obj.kine_states.pos_vecs = out_pos;
            obj.rod.pos_vects = out_pos;
            obj.kine_states.dir_vecs = out_dir;
            obj.rod.dir_vects = out_dir;
        end

        function obj = dynamic_step(obj, dt)
            temp = obj.dynamic_rates(dt);
            obj.dyna_states.rate_vecs = update_dynamics(obj.dyna_states.rate_vecs, temp);
            obj.dyna_states.velocity_vecs = obj.dyna_states.rate_vecs(1,:,:);
            obj.dyna_states.omega_vecs = obj.dyna_states.rate_vecs(2,:,:);
            obj.rod.vel_vects = obj.dyna_states.velocity_vecs;
            obj.rod.omega_vects = obj.dyna_states.omega_vecs;
        end

        function time = sim_step(obj, time, dt)
            obj.kinematic_step(time, dt);
            time = time + obj.prefac(dt);
            obj.update_internal_forces_and_torques();
            obj.dynamic_step(dt);
            obj.kinematic_step(time, dt);
            time = time + obj.prefac(dt);
        end
    end
end

function [out_pos, out_dir] = update_kinematics(prefac, pos_vecs,...
    dir_vecs, velocity_vecs, omega_vecs)
    % updating kinematics
    out_pos = pos_vecs + (prefac .* velocity_vecs);
    rotation_matrix = rodrigues_roth(prefac * omega_vecs, 1.0);
    out_dir = batchMatMat(rotation_matrix, dir_vecs);
end

function rate_vecs = update_dynamics(rate_vecs, scaled_second_deriv_array)
    % updating dynamic_rates
    rate_truncated = rate_vecs(:, :,1:size(scaled_second_deriv_array,3));
    rate_vecs(1:2, :, 1:size(scaled_second_deriv_array,3)) = rate_truncated(1:2, :, 1:size(scaled_second_deriv_array,3)) + scaled_second_deriv_array;
end
