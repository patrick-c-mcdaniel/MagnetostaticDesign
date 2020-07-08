function [ GLz ] = genGLz( L )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
    GLz = zeros(2*L+3, 2*L+1);
    
   
    
    Nr = 2*L+3;
    Nc = 2*L+1;
     k1 = 1:(2*L+1);
    for ik1 = 1:numel(k1)
        k = k1(ik1);
        GLz( k+1, k ) =   sqrt( k*(2*L+2-k) ) / sqrt((2*L+1)*(2*L+3));
    end
    
    
end

