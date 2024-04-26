classdef TargetPoint
    %TARGETPOINT This class carries the location, direction and
    % costs associated with a target

    properties
        cont_pose
        cont_dir
        dis_pose
        dis_dir
        cont_grad_pose
        cont_grad_dir
        dis_grad_pose
        dis_grad_dir
        cost_weight
        target_cost_weight
        pos
        dir
    end

    methods
        function obj = TargetPoint(pos, dir, cost_weight, target_cost_weight)
            obj.cont_pose = 0;
            obj.cont_dir = 0;
            obj.dis_pose = 0;
            obj.dis_dir = 0;
            obj.cont_grad_pose = 0;
            obj.cont_grad_dir = 0;
            obj.dis_grad_pose = zeros(3,1);
            obj.dis_grad_dir = zeros(3,1);
            obj.cost_weight = cost_weight;
            obj.target_cost_weight = target_cost_weight;
            obj.pos = pos;
            obj.dir = dir;
        end

        function obj = update_costs(obj, position, director)
            obj = obj.calculate_discrete_cost_gradient_wrt_position(position);
            obj = obj.calculate_discrete_cost_gradient_wrt_director(director);
        end

        function obj = calculate_discrete_cost_gradient_wrt_position(obj, position)
            % calculate_discrete_cost_gradient_wrt_position.
            % Calculate the gradient of the cost function with respect to the position
            pos_temp = 0.5 * (position(:, end) + position(:, end-1));
            obj.dis_grad_pose(:, end) = ...
                obj.target_cost_weight.position .* (pos_temp - obj.pos);
        end

        function obj = calculate_discrete_cost_gradient_wrt_director(obj, director)
            % calculate_discrete_cost_gradient_wrt_director.
            % Calculate the gradient of the cost function with respect to the director
                dir_temp = director(:, :, end);
                v = zeros(3, 1);
                skew_symmetric_matrix = dir_temp * obj.dir' - obj.dir * dir_temp';
                v(1) = skew_symmetric_matrix(2, 3);
                v(2) = -skew_symmetric_matrix(1, 3);
                v(3) = skew_symmetric_matrix(1, 2);
                obj.dis_grad_dir(:, end) = ...
                obj.target_cost_weight.director .* (dir_temp' * v);
        end

    end
end