classdef kine_state
    %KINE_STATE Summary of this class goes here
    %   Detailed explanation goes here

    properties
        pos_vecs
        dir_vecs
    end

    methods
        function obj = kine_state(pos_vecs,dir_vecs)
        %KINE_STATE Construct an instance of this class
        %   Detailed explanation goes here
        obj.pos_vecs = pos_vecs;
        obj.dir_vecs = dir_vecs;
        end
    end
end