classdef longitudinal_muscle
    %LONGITUDINAL_MUSCLE Implementation of a transverse muscle
    % actuator class

    properties
        nElems
        external_forces
        internal_forces
        internal_couples
        external_couples
        s
        normalized_lengths
        rest_lengths
        lengths
        tangents
        strains
        pos
        ratio_pos
        rest_areas
        areas
        activation
        s_activation
        max_stress
        muscle_force 
        s_force
    end

    methods
        function obj = longitudinal_muscle(nElems, muscle_init_angle, ...
                ratio_position, rest_area, max_stress)

            obj.nElems = nElems;
            transform_matrix = [cos(muscle_init_angle); sin(muscle_init_angle); 0];
            obj.ratio_pos =  ratio_position .* transform_matrix;
            obj.s = linspace(0, 1, nElems + 1);
            obj.normalized_lengths = zeros(1,nElems);
            obj.rest_lengths = zeros(1,nElems);
            obj.lengths = zeros(1,nElems);
            obj.tangents = zeros(3,nElems);
            obj.strains = zeros(3,nElems);
            obj.pos = zeros(3,nElems);

            obj.rest_areas = rest_area;
            obj.areas = obj.rest_areas;
            obj.internal_forces = zeros(3, nElems);  % material frame
            obj.external_forces = zeros(3, nElems + 1);  % global frame
            obj.internal_couples = zeros(3, nElems - 1);  % material frame
            obj.external_couples = zeros(3, nElems);  % material frame
            obj.activation = zeros(1,nElems);
            obj.s_activation = (obj.s(1:end-1) + obj.s(2:end)) / 2;
            obj.max_stress = max_stress;
            obj.muscle_force = zeros(1, obj.nElems);
            obj.s_force = 0.5 * (obj.s(1:end-1) + obj.s(2:end));
        end

        function obj = update_muscle_kyn_shear(obj, rod)

            % Calc muscle area
            obj.areas = obj.rest_areas ./ (rod.dilatations + 1e-14);
            % Calc muscle pos
            obj.pos = rod.radius .* obj.ratio_pos;
            % Calc muscle strain
            % Initialize the shear matrix with the same size as sigma
            shear = rod.sigma;
            shear(3, :) = rod.sigma(3, :) + 1;
            pos_derivative = diff(obj.pos, 1, 2) ./ (rod.rest_voronoi_lengths .* rod.voronoi_dilatations + 1e-14);

            obj.strains = shear + trapezoidal_rule(cross(rod.kappa, average2D(obj.pos))...
                            + pos_derivative);
            % Calc muscle tangents
            magnitudes = sqrt(sum(obj.strains .^ 2, 1));
            
            obj.tangents = obj.strains ./ (magnitudes + + 1e-14);
        end

        function obj = update_muscle_force_couple(obj, rod)
            obj = obj.update_muscle_kyn_shear(rod);
            obj.lengths = sqrt(sum(obj.strains .^ 2, 1));
            obj = obj.calc_norm_length();

            % Calc muscle force.
            obj.muscle_force = (obj.activation .* obj.max_stress .* force_length_weight(obj.normalized_lengths))...
                .* obj.areas;
            
            % Calc internal and external forces / couples.
            obj = obj.calc_force_couple(rod);

        end

        function obj = set_length_as_rest(obj, rod)
            obj = obj.update_muscle_kyn_shear(rod);
            obj.lengths = sqrt(sum(obj.strains .^ 2, 1));
            obj.rest_lengths = obj.lengths;
        end

        function obj = calc_norm_length(obj)
            obj.normalized_lengths = obj.lengths ./ (obj.rest_lengths + 1e-14);
        end

        function obj = calc_force_couple(obj, rod)
            obj.internal_forces = obj.muscle_force .* obj.tangents;
            obj.internal_couples = average2D(cross(obj.pos, obj.internal_forces));
            obj = obj.internal2externalLoad(rod);
        end

        function obj = reset_actuations(obj)
            obj.internal_couples = obj.internal_couples .* 0;
            obj.external_couples = obj.external_couples .* 0;
            obj.internal_forces = obj.internal_forces .* 0;
            obj.external_forces = obj.external_forces .* 0;
        end

        function obj = set_activation(obj, activations)
            obj.activation = activations;
        end

        function obj = internal2externalLoad(obj, cosserat_rod)

            obj.external_forces = two_point_rule(...
                material2labframe(cosserat_rod.dir_vects, obj.internal_forces),[]);
        
            obj.external_couples = (two_point_rule(obj.internal_couples) +...
                    trapezoidal_rule( cross(cosserat_rod.kappa, obj.internal_couples,1) ...
                    .* cosserat_rod.rest_voronoi_lengths) + cross(...
                        batchMatVec(cosserat_rod.dir_vects, cosserat_rod.tangents .* cosserat_rod.dilatations),...
                        obj.internal_forces,1) .* cosserat_rod.rest_lengths);
        end
    end
end