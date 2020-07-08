function [ sz3tps ] = Mxyz_to_sz3tps( Mxyz, sss )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    Mx = Mxyz(1:3:end);
    My = Mxyz(2:3:end);
    Mz = Mxyz(3:3:end);
    
    Mmax = 0.0254^3*1.42/(4*pi*1e-7);
    
    Mm  = sqrt( Mx.^2 + My.^2 + Mz.^2 );
    
    tts = acos(Mz./Mm);
    tts(Mm==0) = 0;
    
    pps = atan2(My, Mx);
    tps = [tts(:) pps(:) sss(:)];
    
    Sx = 0.0254*ones(size(Mx));
    Sz = 0.0254*ones(size(Mx));
    Sy = Mm./Mmax*0.0254;
    
    sz3 = [Sx(:) Sy(:) Sz(:)];
    
    sz3tps = [sz3 tps];
end

