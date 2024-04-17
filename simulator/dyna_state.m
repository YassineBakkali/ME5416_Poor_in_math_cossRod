classdef dyna_state
    %DYNA_STATE Summary of this class goes here
    %   Detailed explanation goes here

    properties
        rate_vecs
        dvdt_dwdt_vecs
        velocity_vecs
        omega_vecs
    end

    methods
        function obj = dyna_state(rate_vecs,dvdt_dwdt_vecs,velocity_vecs,...
                omega_vecs)
        %DYNA_STATE Construct an instance of this class
        %   Detailed explanation goes here
        obj.rate_vecs = rate_vecs;
        obj.dvdt_dwdt_vecs = dvdt_dwdt_vecs;
        obj.velocity_vecs = velocity_vecs;
        obj.omega_vecs = omega_vecs;
        end

        function [velocity_rates, omega_rates] = kinematic_rates(obj)
        %METHOD1 Summary of this method goes here
        %   Detailed explanation goes here
        velocity_rates = obj.velocity_vecs;
        omega_rates = obj.omega_vecs;
        end

        function [acc_and_alpha] = dynamic_rates(obj, prefac)
        %METHOD1 Summary of this method goes here
        %   Detailed explanation goes here
        acc_and_alpha = prefac .* obj.dvdt_dwdt_vecs;
        end
    end
end