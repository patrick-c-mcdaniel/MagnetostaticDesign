function [ bx by bz ] = vec2field( bvec )
% vec2field - turns Dipole matrix output vector of B-field components into
% Bx, By, and Bz component vectors (accounts for structure of dipole
% matrix)

    bx = bvec(1:3:end);
    by = bvec(2:3:end);
    bz = bvec(3:3:end);


end

