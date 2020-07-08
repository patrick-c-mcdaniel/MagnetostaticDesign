function [ GLy ] = genGLy( L )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
    GLy = zeros(2*L+3, 2*L+1);
    
   
    
    Nr = 2*L+3;
    Nc = 2*L+1;
     k1 = 1:(L-1);
    for ik1 = 1:numel(k1)
        k = k1(ik1);
        GLy( 2*L+2-k, k       ) =   sqrt( k*(k+1))/2/sqrt((2*L+1)*(2*L+3));
        GLy( k+2,     2*L+2-k ) = -1*sqrt( k*(k+1))/2/sqrt((2*L+1)*(2*L+3));      
    end
    
    k2 = 1:(L);
    for ik2 = 1:numel(k2)
        k = k2(ik2);
        GLy( k,       2*L+2-k ) =  -1*sqrt((2*L+2-k)*(2*L+3-k))/2/sqrt((2*L+1)*(2*L+3));
        GLy( 2*L+4-k, k       ) =     sqrt((2*L+2-k)*(2*L+3-k))/2/sqrt((2*L+1)*(2*L+3));
    end
    
    GLy(L+1, L+1) = -1*sqrt( 2*(L+1)*(L+2) ) / 2/sqrt( (2*L+1)*(2*L+3) );
    GLy(L+2, L  ) = sqrt( 2*L*(L+1))/ 2/sqrt( (2*L+1)*(2*L+3) );
    
end

