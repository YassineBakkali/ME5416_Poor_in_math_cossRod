function [ss_mat] = skew_sym(in)
%SKEW_SYM Summary of this function goes here
    % This function takes a 3x1 vector and returns a skew-symmetric matrix
    if length(in) ~= 3
        error('Input must be a 3-element vector.');
    end
    
    ss_mat = [  0   -in(3)   in(2);
          in(3)   0   -in(1);
         -in(2)  in(1)    0];
end

