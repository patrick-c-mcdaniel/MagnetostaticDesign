function [ Bxyz ] = Bcomp135_sz3_tps_1var_debug( D135, J3, J5, sz3tps, errcoefs )
% Compute B-field from sizes+orientations using numerically-derived block SH error coeffs
%   - Converts Mxyz (dipole terms) to (Mxyz + M5) (dipole + L=5 multipole)
%   - Calculates B-field using L=1 + L=5 terms
%   - Assumes block dipole moment varies exclusively by changing
%       y-dimension (ie the x- and z- dimensions of all magnet blocks are
%       fixed, generally at 0.0254m. 
%
%   Inputs:
%       D135     : Field comp. matrix (L=1,3,5 matrices concatenated
%                   together)
%       J3       : Permutation matrix for L=3 SH terms corresponding to
%                   swapping the y and z axes (for fast calculation of L=3 unitary
%                   transformation)
%       J5       : Permutation matrix for L=5 SH terms corresponding to
%                   swapping the y and z axes (for fast calculation of L=5 unitary
%                   transformation)
%       sz3      : (Nbk)x(3) matrix of (x)x(y)x(z) dimensions of blocks
%       tps      : (Nbk)x(3) matrix of (theta, phi, psi) orientations of
%                   all blocks
%       effcoefs : error coeffs for (l;m)=[ 1 3 3 5 5 5; 0 0 2 0 2 4 ].
%                   Should be an (Nterm)x(Npoly) = (6)x(6) matrix
%
%   Outputs:
%       Bxyz    : B-field at all target points


    %% convert Mxyz (dipole components) to M5 (L=5 multipole components)
    sz3 = sz3tps(:,1:3);
    tps = sz3tps(:,4:6);
    
    %%% scale down to be realistic...
%     Br_des = 1.6593;
    Br = 1.42;

    tts = tps(:,1);
    pps = tps(:,2);
    sss = tps(:,3);
    
    
    MUZ = 4*pi*1e-7;
    Mblk = Br/MUZ;
    
    Nbk = size(sz3,1);
    
    Xs1 = sz3(:,1);
    Xs3 = Xs1.^3;
    Xs5 = Xs1.^5;
    
    Zs1 = sz3(:,3);
    Zs3 = Zs1.^3;
    Zs5 = Zs1.^5;
    
    Ys1 = sz3(:,2);
    Ys3 = Ys1.^3;
    Ys5 = Ys1.^5;
    Ys50 = Ys1.^(5:-1:0);
    
    Mm = Mblk * Xs1 .* Ys1 .* Zs1;
    Mx = cos(pps).*sin(tts).*Mm;
    My = sin(pps).*sin(tts).*Mm;
    Mz = cos(tts).*Mm;
    Mxyz = [ Mx(:), My(:), Mz(:)]';
    Mxyz = Mxyz(:);
    
    %% manage error terms
    M1off = zeros(Nbk,1);
    M3err = zeros(Nbk*7,1);
    M5err = zeros(Nbk*11,1);
    
    M1off(1+1*(0:295))   = sum(repmat(errcoefs(1,:),[Nbk, 1]) .* Ys50, 2);
    M3err(4+7*(0:295))   = sum(repmat(errcoefs(2,:),[Nbk, 1]) .* Ys50, 2);
    M3err(6+7*(0:295))   = sum(repmat(errcoefs(3,:),[Nbk, 1]) .* Ys50, 2);
    M5err(6+11*(0:295))  = sum(repmat(errcoefs(4,:),[Nbk, 1]) .* Ys50, 2);
    M5err(8+11*(0:295))  = sum(repmat(errcoefs(5,:),[Nbk, 1]) .* Ys50, 2);
    M5err(10+11*(0:295)) = sum(repmat(errcoefs(6,:),[Nbk, 1]) .* Ys50, 2);
    
    %% analytical multipole terms
    
    M3     = zeros(numel(Mxyz)*7/3, 1);
    M3rot  = zeros(numel(Mxyz)*7/3, 1);
    M5     = zeros(numel(Mxyz)*11/3, 1);
    M5rot  = zeros(numel(Mxyz)*11/3, 1);
    
    %%% analytical
    M3(4+7*(0:295)) = 1/4*sqrt(7/pi)   * Mblk * (1/2*Xs1.*Ys1.*Zs3 - 1/4*Xs3.*Ys1.*Zs1 - 1/4*Xs1.*Ys3.*Zs1);
    M3(6+7*(0:295)) = 1/4*sqrt(105/pi) * Mblk * (1/12*Xs3.*Ys1.*Zs1 - 1/12*Xs1.*Ys3.*Zs1);
    
    M5(6+11*(0:295))  = 1/16*sqrt(11/pi)  *Mblk*(1/2*Xs1.*Ys1.*Zs5  - 5/6*Xs3.*Ys1.*Zs3  - 5/6*Xs1.*Ys3.*Zs3  + 3/16*Xs5.*Ys1.*Zs1 + 5/24*Xs3.*Ys3.*Zs1 + 3/16*Xs1.*Ys5.*Zs1);
    M5(8+11*(0:295))  = 1/8 *sqrt(1155/pi)*Mblk*(1/24*Xs3.*Ys1.*Zs3 - 1/24*Xs1.*Ys3.*Zs3 - 1/80*Xs5.*Ys1.*Zs1 + 1/80*Xs1.*Ys5.*Zs1);
    M5(10+11*(0:295)) = 3/16*sqrt(385/pi) *Mblk*(1/80*Xs5.*Ys1.*Zs1 + 1/80*Xs1.*Ys5.*Zs1 - 1/24*Xs3.*Ys3.*Zs1);
    
    
    M1sc = repmat( (Mm+M1off)./Mm, [1 3]);
    M1sc = M1sc';
    M1sc = M1sc(:);
    
    M1 = Mxyz.*M1sc;
    M3 = M3 + M3err;
    M5 = M5 + M5err;
    
% pause
    for ibb = 1:numel(Mxyz)/3
        pp = pps(ibb);
        tt = tts(ibb);
        ss = sss(ibb);

        Mtmp3 = M3((1+(ibb-1)*7):(ibb*7));
        M3rot((1+(ibb-1)*7):(ibb*7)) = Xsh3(pp)*J3*Xsh3(tt)*J3*Xsh3(ss)*Mtmp3(:);
        
        Mtmp5 = M5((1+(ibb-1)*11):(ibb*11));
        M5rot((1+(ibb-1)*11):(ibb*11)) = Xsh5(pp)*J5*Xsh5(tt)*J5*Xsh5(ss)*Mtmp5(:);

    end

    
    %% compute B-field
    
    Bxyz =  D135 * [ M1; M3rot; M5rot ];
        
    M135 = [ M1; M3rot; M5rot ];
    save('M135_full_debut.mat','M135');
%     Bxyz =  D135 * [ M1; zeros(size(M3rot)); zeros(size(M5rot)) ];
    
   
end

