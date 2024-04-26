classdef kine_state


    properties
        pos_vecs
        dir_vecs
    end

    methods
        function obj = kine_state(pos_vecs,dir_vecs)

        obj.pos_vecs = pos_vecs;
        obj.dir_vecs = dir_vecs;
        end
    end
end