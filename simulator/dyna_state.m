classdef dyna_state

    properties
        rate_vecs
        dvdt_dwdt_vecs
        velocity_vecs
        omega_vecs
    end

    methods
        function obj = dyna_state(rate_vecs,dvdt_dwdt_vecs,velocity_vecs,...
                omega_vecs)

        obj.rate_vecs = rate_vecs;
        obj.dvdt_dwdt_vecs = dvdt_dwdt_vecs;
        obj.velocity_vecs = velocity_vecs;
        obj.omega_vecs = omega_vecs;
        end

        function [velocity_rates, omega_rates] = kinematic_rates(obj)

        velocity_rates = obj.velocity_vecs;
        omega_rates = obj.omega_vecs;
        end

        function [acc_and_alpha] = dynamic_rates(obj, prefac)

        acc_and_alpha = prefac .* obj.dvdt_dwdt_vecs;
        end
    end
end