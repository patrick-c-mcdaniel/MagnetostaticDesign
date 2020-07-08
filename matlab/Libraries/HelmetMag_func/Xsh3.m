function [ X ] = Xsh3( ANG )
%XSH Summary of this function goes here
%   Detailed explanation goes here
    X =      diag([ cos(3*ANG),  cos(2*ANG),  cos(ANG),  1, cos(ANG), cos(2*ANG), cos(3*ANG)]) + ...
        flip(diag([ sin(-3*ANG), sin(-2*ANG), sin(-ANG), 0, sin(ANG), sin(2*ANG), sin(3*ANG)]), 1); 

end

