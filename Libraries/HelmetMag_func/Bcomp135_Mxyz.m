function [ Bxyz ] = Bcomp135_Mxyz( D135, J3, J5, Mred, sssred, errcoefs, Bxyz_off )
% Compute B-field from Mred + psi's using numerically-derived block SH 
%   error coeffs; exploit symmetry of magnet for faster computation
%   - Converts Mred (reduced set dipoles) to reduced multipole set
%   - Compute B-field using reduced D135 matrix
%   - Assumes block dipole moment varies exclusively by changing
%       y-dimension (ie the x- and z- dimensions of all magnet blocks are
%       fixed, generally at 0.0254m. 
%
%   Inputs:
%       D135red  : Field comp. matrix (L=1,3,5 matrices concatenated
%                   together; compressed to exploit symmetry)
%       J3       : Permutation matrix for L=3 SH terms corresponding to
%                   swapping the y and z axes (for fast calculation of L=3 unitary
%                   transformation)
%       J5       : Permutation matrix for L=5 SH terms corresponding to
%                   swapping the y and z axes (for fast calculation of L=5 unitary
%                   transformation)
%       Mred     : (Nred)x(1) vector of reduced set dipole components
%       sssred   : reduced set of block psi angles
%       Nmid     : Number of midline blocks (listed first in Mred - only 2
%                   components since My=0 along midline)
%       effcoefs : error coeffs for (l;m)=[ 1 3 3 5 5 5; 0 0 2 0 2 4 ].
%                   Should be an (Nterm)x(Npoly) = (6)x(6) matrix
%       Bxyz_off : offset field (ie a background field term due to coil;
%                   other magnet; etc.). Vector with the same format as
%                   Bxyz
%   Outputs:
%       Bxyz    : B-field at all target points


    %% convert Mxyz (dipole components) to M5 (L=5 multipole components)
    
    sz3tps = Mxyz_to_sz3tps( Mred, sssred );
    
    Nbr = numel(sssred);
    
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
    Ys50(Ys1==0,end) = 0;
    
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
    
    M1off(1+1*(0:(Nbr-1)))   = sum(repmat(errcoefs(1,:),[Nbk, 1]) .* Ys50, 2);
    M3err(4+7*(0:(Nbr-1)))   = sum(repmat(errcoefs(2,:),[Nbk, 1]) .* Ys50, 2);
    M3err(6+7*(0:(Nbr-1)))   = sum(repmat(errcoefs(3,:),[Nbk, 1]) .* Ys50, 2);
    M5err(6+11*(0:(Nbr-1)))  = sum(repmat(errcoefs(4,:),[Nbk, 1]) .* Ys50, 2);
    M5err(8+11*(0:(Nbr-1)))  = sum(repmat(errcoefs(5,:),[Nbk, 1]) .* Ys50, 2);
    M5err(10+11*(0:(Nbr-1))) = sum(repmat(errcoefs(6,:),[Nbk, 1]) .* Ys50, 2);
    
    %% analytical multipole terms
    
    M3     = zeros(numel(Mxyz)*7/3, 1);
    M3rot  = zeros(numel(Mxyz)*7/3, 1);
    M5     = zeros(numel(Mxyz)*11/3, 1);
    M5rot  = zeros(numel(Mxyz)*11/3, 1);
    
    %%% analytical
    M3(4+7*(0:(Nbr-1))) = 1/4*sqrt(7/pi)   * Mblk * (1/2*Xs1.*Ys1.*Zs3 - 1/4*Xs3.*Ys1.*Zs1 - 1/4*Xs1.*Ys3.*Zs1);
    M3(6+7*(0:(Nbr-1))) = 1/4*sqrt(105/pi) * Mblk * (1/12*Xs3.*Ys1.*Zs1 - 1/12*Xs1.*Ys3.*Zs1);
    
    M5(6+11*(0:(Nbr-1)))  = 1/16*sqrt(11/pi)  *Mblk*(1/2*Xs1.*Ys1.*Zs5  - 5/6*Xs3.*Ys1.*Zs3  - 5/6*Xs1.*Ys3.*Zs3  + 3/16*Xs5.*Ys1.*Zs1 + 5/24*Xs3.*Ys3.*Zs1 + 3/16*Xs1.*Ys5.*Zs1);
    M5(8+11*(0:(Nbr-1)))  = 1/8 *sqrt(1155/pi)*Mblk*(1/24*Xs3.*Ys1.*Zs3 - 1/24*Xs1.*Ys3.*Zs3 - 1/80*Xs5.*Ys1.*Zs1 + 1/80*Xs1.*Ys5.*Zs1);
    M5(10+11*(0:(Nbr-1))) = 3/16*sqrt(385/pi) *Mblk*(1/80*Xs5.*Ys1.*Zs1 + 1/80*Xs1.*Ys5.*Zs1 - 1/24*Xs3.*Ys3.*Zs1);
    
    M1sc_1d = (Mm+M1off)./Mm;
    M1sc_1d(Mm==0) = 0;
    
    M1sc = repmat( M1sc_1d, [1 3]);
%     M1sc(Mm==0) = 0;
    
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

%     M135_red = [ M1; M3rot; M5rot ];
    
    
    %% compute B-field
   
    Bxyz =  D135 * [ M1; M3rot; M5rot ];
    
    if nargin<7
        Bxyz_off = zeros(size(Bxyz));
    end
    Bxyz = Bxyz + Bxyz_off;
    
   
end

