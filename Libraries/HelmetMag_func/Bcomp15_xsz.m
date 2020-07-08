function [ Bxyz ] = Bcomp135_xsz( D135, J5, Mxyz )
% Compute B-field from Mxyz
%   - Converts Mxyz (dipole terms) to (Mxyz + M5) (dipole + L=5 multipole)
%   - Calculates B-field using L=1 + L=5 terms
%   - Assumes block dipole moment varies exclusively by changing
%       x-dimension (ie the y- and z- dimensions of all magnet blocks are
%       fixed, generally at 0.0254m. 
%
%   Inputs:
%       D15     : Field comp. matrix (L=1 + L=5 matrices concatenated
%                   together)
%       J5      : Permutation matrix for L=5 SH terms corresponding to
%                   swapping the y and z axes (for fast calculation of L=5 unitary
%                   transformation)
%       Mxyz    : Dipole moment vector (ie description of magnet design)
%
%   Outputs:
%       Bxyz    : B-field at all target points


    %% convert Mxyz (dipole components) to M5 (L=5 multipole components)
    
    Mx = Mxyz(1:3:end);
    My = Mxyz(2:3:end);
    Mz = Mxyz(3:3:end);

    mms = sqrt( Mx.^2 + My.^2 + Mz.^2 );
    tts = acos( Mz ./ mms );
    pps = atan2( My, Mx );

    Br = 1.6593;
    MUZ = 4*pi*1e-7;
    Mblk = Br/MUZ;
    
    Ysz = 0.0254;
    Zsz = 0.0254;
    Xsz = mms / Mblk / Ysz / Zsz;
    Xs3 = Xsz.^3;
    Xs5 = Xsz.^5;
    
    M3     = zeros(numel(Mxyz)*7/3, 1);
    M3rot  = zeros(numel(Mxyz)*7/3, 1);
    M5     = zeros(numel(Mxyz)*11/3, 1);
    M5rot  = zeros(numel(Mxyz)*11/3, 1);
    
    M3(4+7*(0:295)) = 1/4*sqrt(7/pi) * Mblk * (1/2*Xsz*Ysz*Zsz^3 - 1/4*Xs3*Ysz*Zsz - 1/4*Xsz*Ysz^3*Zsz);



    %%% analytical
    M5(6+11*(0:295)) = (-0.0728 + (0.0005)/2) * m5s;
    % M5(2+11*(0:295))  = -0.3 * m5s;
    M5(10+11*(0:295)) = (-0.0367 + (-0.0005/2)) * m5s;

    for ibb = 1:numel(Mxyz)/3
        pp = pps(ibb);
        tt = tts(ibb);

        Mtmp5 = M5((1+(ibb-1)*11):(ibb*11));
        M5rot((1+(ibb-1)*11):(ibb*11)) = Xsh5(pp)*J5*Xsh5(tt)*J5*Mtmp5(:);

    end

    
    %% compute B-field
    
    Bxyz =  D15 * [ Mxyz; M5rot ];
    
   
end

