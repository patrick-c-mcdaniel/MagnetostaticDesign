function [ sz3tps ] = Mred_to_sz3tps( Mred, sssred, Nmid )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    Mxmid = Mred(1:2:(2*Nmid));
    Mymid = zeros(size(Mxmid));
    Mzmid = Mred(2:2:(2*Nmid));
    
    Mxpos = Mred( (2*Nmid+1):3:end );
    Mypos = Mred( (2*Nmid+2):3:end );
    Mzpos = Mred( (2*Nmid+3):3:end );
    
    Mx = cat(1, Mxmid(:), Mxpos(:));
    My = cat(1, Mymid(:), Mypos(:));
    Mz = cat(1, Mzmid(:), Mzpos(:));
    
    Mmax = 0.0254^3*1.42/(4*pi*1e-7);
    
  
    Mm  = sqrt( Mx.^2 + My.^2 + Mz.^2 );
    
    tts = acos(Mz./Mm);
    tts(Mm==0) = 0;
    
    pps = atan2(My, Mx);
    
    tps = [tts(:) pps(:) sssred(:)];
    
    Sx = 0.0254*ones(size(Mx));
    Sz = 0.0254*ones(size(Mx));
    Sy = Mm./Mmax*0.0254;
    
    sz3 = [Sx(:) Sy(:) Sz(:)];
   
    sz3tps = [sz3 tps];
end

