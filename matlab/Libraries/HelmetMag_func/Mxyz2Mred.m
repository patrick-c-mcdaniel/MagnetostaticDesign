function [ Mred ] = Mxyz2Mred( Mxyz, Imid, Iypos )
% converts full Mxyz vector to reduced Mred vector. Mred is a compressed,
% symmetrized description of Mxyz; use of Mred imposes symmetry on the
% magnet design and reduces the size of the field computation/optimization
% problem
%
% Inputs:
%       Mxyz    - Full magnet vector; size (3*Nmag)x(1). Entries 1:3:end ->
%                   x-components, etc.
%       Imid    - Indices of magnets (pptx convention) located along the
%                   XZ-plane (symmetry plane). Keeps x- and z- components
%                   for these blocks; y-components forced to be ==0;
%       Iypos   - Indices of blocks with positive y-coordinate (must be
%                   sorted in ascending order). 
%       Iyneg   - Indices of blocks with negative y-coordinate (must
%                   pairwise match Iypos)
%
% Outputs:
%       Mred    - Reduced magnet vector; order goes: { Mmid_x1, Mmid_z1,
%                   Mmid_x2, Mmid_z2, Mmid_x3, ... Mmid_zNmid, Mpos_x1, Mpos_y1,
%                   Mpos_z1, Mpos_x2, ... Mpos_zNpos }

    Nmid = numel(Imid);
    Npos = numel(Iypos);
    
    Mred = zeros( 2*Nmid + 3*Npos, 1 );
    Mred(1:2:(2*Nmid)) = Mxyz(3*Imid-2);
    Mred(2:2:(2*Nmid)) = Mxyz(3*Imid-0);
    
    Mred(2*Nmid + (1:3:(3*Npos))) = Mxyz(3*Iypos-2);
    Mred(2*Nmid + (2:3:(3*Npos))) = Mxyz(3*Iypos-1);
    Mred(2*Nmid + (3:3:(3*Npos))) = Mxyz(3*Iypos-0);

end

