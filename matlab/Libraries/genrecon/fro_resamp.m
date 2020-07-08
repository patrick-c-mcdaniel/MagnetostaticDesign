function [ fro_resamp ] = fro_resamp( a_struct )
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here
    
    phant = a_struct.phant_struct;
    xphant = phant.xpts;
    
    fro_resamp = interp1( a_struct.ro_struct.xpos, a_struct.ro_struct.freqs, xphant, 'cubic' );

end

