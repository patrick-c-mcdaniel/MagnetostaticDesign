function [ Mxyz ] = sz3tps_to_Mxyz( sz3tps )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    Sx = sz3tps(:,1);
    Sz = sz3tps(:,3);
    Sy = sz3tps(:,2);
    
    tts = sz3tps(:,4);
    pps = sz3tps(:,5);
    sss = sz3tps(:,6);
    
    M = 1.42/(4*pi*1e-7);
    Mm = M.*Sx.*Sy.*Sz;

    Mx = Mm.*sin(tts).*cos(pps-pi/2);
    My = Mm.*sin(tts).*sin(pps-pi/2);
    Mz = Mm.*cos(tts);
    
    Mxyz = [Mx(:) My(:) Mz(:)];
    Mxyz = Mxyz';
    Mxyz = Mxyz(:);
    

end

