function [ Bxyz ] = Bcomp135_xsz( D135, J3, J5, Mxyz )
% Compute B-field from Mxyz
%   - Converts Mxyz (dipole terms) to (Mxyz + M5) (dipole + L=5 multipole)
%   - Calculates B-field using L=1 + L=5 terms
%   - Assumes block dipole moment varies exclusively by changing
%       y-dimension (ie the x- and z- dimensions of all magnet blocks are
%       fixed, generally at 0.0254m. 
%
%   Inputs:
%       D135    : Field comp. matrix (L=1,3,5 matrices concatenated
%                   together)
%       J3      : Permutation matrix for L=3 SH terms corresponding to
%                   swapping the y and z axes (for fast calculation of L=3 unitary
%                   transformation)
%       J5      : Permutation matrix for L=5 SH terms corresponding to
%                   swapping the y and z axes (for fast calculation of L=5 unitary
%                   transformation)
%       Mxyz    : Dipole moment vector (ie description of magnet design)
%
%   Outputs:
%       Bxyz    : B-field at all target points


    %% convert Mxyz (dipole components) to M5 (L=5 multipole components)
    
    
    %%% scale down to be realistic...
    Br_des = 1.6593;
    Br = 1.42;
    Mxyz = Mxyz * Br/Br_des;
    
%     %%% unrealistic (design) Br
%     Br = 1.6593;
    
    
    Mx = Mxyz(1:3:end);
    My = Mxyz(2:3:end);
    Mz = Mxyz(3:3:end);

    mms = sqrt( Mx.^2 + My.^2 + Mz.^2 );
    tts = acos( Mz ./ mms );
    pps = atan2( My, Mx );


    
    
    
    MUZ = 4*pi*1e-7;
    Mblk = Br/MUZ;
    
    Xs1 = 0.0254;
    Xs3 = Xs1^3;
    Xs5 = Xs1^5;
    
    Zs1 = 0.0254;
    Zs3 = Zs1^3;
    Zs5 = Zs1^5;
    
    Ys1 = mms / Mblk / Xs1 / Zs1;
    Ys3 = Ys1.^3;
    Ys5 = Ys1.^5;
    
    M3     = zeros(numel(Mxyz)*7/3, 1);
    M3rot  = zeros(numel(Mxyz)*7/3, 1);
    M5     = zeros(numel(Mxyz)*11/3, 1);
    M5rot  = zeros(numel(Mxyz)*11/3, 1);
    
    %%% analytical
    M3(4+7*(0:295)) = 1/4*sqrt(7/pi)   * Mblk * (1/2*Xs1*Ys1*Zs3 - 1/4*Xs3*Ys1*Zs1 - 1/4*Xs1*Ys3*Zs1);
    M3(6+7*(0:295)) = 1/4*sqrt(105/pi) * Mblk * (1/6*Xs3*Ys1*Zs1 - 1/6*Xs1*Ys3*Zs1);
    
    M5(6+11*(0:295))  = 1/16*sqrt(11/pi)  *Mblk*(1/2*Xs1*Ys1*Zs5  - 5/6*Xs3*Ys1*Zs3  - 5/6*Xs1*Ys3*Zs3  + 3/16*Xs5*Ys1*Zs1 + 5/24*Xs3*Ys3*Zs1 + 3/16*Xs1*Ys5*Zs1);
    M5(8+11*(0:295))  = 1/8 *sqrt(1155/pi)*Mblk*(1/24*Xs3*Ys1*Zs3 - 1/24*Xs1*Ys3*Zs3 - 1/80*Xs5*Ys1*Zs1 + 1/80*Xs1*Ys5*Zs1);
    M5(10+11*(0:295)) = 3/16*sqrt(385/pi) *Mblk*(1/80*Xs5*Ys1*Zs1 + 1/80*Xs1*Ys5*Zs1 - 1/24*Xs3*Ys3*Zs1);
% pause
    for ibb = 1:numel(Mxyz)/3
        pp = pps(ibb);
        tt = tts(ibb);

        Mtmp3 = M3((1+(ibb-1)*7):(ibb*7));
        M3rot((1+(ibb-1)*7):(ibb*7)) = Xsh3(pp)*J3*Xsh3(tt)*J3*Mtmp3(:);
        
        Mtmp5 = M5((1+(ibb-1)*11):(ibb*11));
        M5rot((1+(ibb-1)*11):(ibb*11)) = Xsh5(pp)*J5*Xsh5(tt)*J5*Mtmp5(:);

    end

    
    %% compute B-field
    
    Bxyz =  D135 * [ Mxyz; M3rot; M5rot ];
    
   
end

