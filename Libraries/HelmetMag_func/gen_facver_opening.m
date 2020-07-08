function [ fac_blk, ver_blk ] = gen_facver_opening( xyz_ctr, xyz_sz, phi_blk, tht_blk, psi_blk, facnum, tcoffer )

%% matrices about connectivity; face numbering, etc.
    fac2ver = [ 1, 2, 3, 4; ...
                3, 4, 7, 8; ...
                5, 6, 7, 8; ...
                1, 2, 5, 6; ...
                1, 4, 5, 8; ...
                2, 3, 6, 7 ];

%%
    xctr = xyz_ctr(1);
    yctr = xyz_ctr(2);
    zctr = xyz_ctr(3);
    
    xsz = xyz_sz(1)/2 + tcoffer;
    ysz = xyz_sz(2)/2 + tcoffer;
    zsz = xyz_sz(3)/2 + tcoffer;
    
    xsz_op = xyz_sz(1)/2;
    ysz_op = xyz_sz(2)/2;
    zsz_op = xyz_sz(3)/2;

    cph = cos(phi_blk);
    cth = cos(tht_blk);
    cps = cos(psi_blk);
    
    sph = sin(phi_blk);
    sth = sin(tht_blk);
    sps = sin(psi_blk);
    
    Rot_eul = [ (cph*cps-cth*sph*sps),      (-1*cph*sps-cth*cps*sph),   sph*sth;        ...
                cps*sph+cph*cth*sps,        (cph*cth*cps-sph*sps),      -1*cph*sth;     ...
                sth*sps,                    cps*sth,                    cth             ];
            
    xver_all = [ -xsz, -xsz, +xsz, +xsz, -xsz, -xsz, +xsz, +xsz ]';
    yver_all = [ -ysz, +ysz, +ysz, -ysz, -ysz, +ysz, +ysz, -ysz ]';
    zver_all = [ +zsz, +zsz, +zsz, +zsz, -zsz, -zsz, -zsz, -zsz ]';
    
    switch facnum
        case 0
        case 1
            zver_all(fac2ver(facnum, :)) = zsz_op;
        case 2
            xver_all(fac2ver(facnum, :)) = xsz_op;
        case 3
            zver_all(fac2ver(facnum, :)) = -zsz_op;
        case 4
            xver_all(fac2ver(facnum, :)) = -xsz_op;
        case 5
            yver_all(fac2ver(facnum, :)) = -ysz_op;
        case 6
            yver_all(fac2ver(facnum, :)) = ysz_op;
        otherwise
            error('invalid face number; must be between 1 and 6');
    end
    
    xyz_ver = cat(2, xver_all, yver_all, zver_all);
    
    ver_corners = Rot_eul * xyz_ver';
    ver_blk = ver_corners' + [ xctr*ones(8,1), yctr*ones(8,1), zctr*ones(8,1)];
    
    fac_blk = [ 1 4 2 ; ...
                3 2 4 ; ...
                4 8 3 ; ...
                7 3 8 ; ...
                8 5 7 ; ...
                6 7 5 ; ...
                5 1 6 ; ...
                2 6 1 ; ...
                5 8 1 ; ...
                4 1 8 ; ...
                6 2 7 ; ...
                3 7 2 ];

end