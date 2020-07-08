function [ ll, pp, rr, RMAX ] = xyz_to_lpr( xx, yy, zz, xyzor )

    ll = zeros(size(xx));
    pp = zeros(size(xx));
    rr = zeros(size(xx));
    
    %%% split into "dome" and "cylinder" points
    
    xx_d = xx(zz<=xyzor(3)) - xyzor(1); xx_c = xx(zz>xyzor(3)) - xyzor(1);
    yy_d = yy(zz<=xyzor(3)) - xyzor(2); yy_c = yy(zz>xyzor(3)) - xyzor(2);
    zz_d = zz(zz<=xyzor(3)) - xyzor(3); zz_c = zz(zz>xyzor(3)) - xyzor(3);
    
    %%% convert "dome" coordinates
    tht_d = atan2( sqrt((xx_d.^2 + yy_d.^2)), abs(zz_d) );
    phi_d = atan2( xx_d, yy_d );
    rad_d = sqrt( xx_d.^2 + yy_d.^2 + zz_d.^2 );
    RMAX = max(rad_d(:));
    
    ll_d = tht_d * RMAX;
    rr_d = rad_d;
%     rr_d = RMAX;
    pp_d = phi_d;
    
    %%% convert "cylinder" coordinates
    ll_c  = pi/2*RMAX + (zz_c);
    phi_c = atan2( xx_c, yy_c );
    rad_c = sqrt( xx_c.^2 + yy_c.^2);
%    rad_c = RMAX;
    
    rr_c = rad_c;
    pp_c = phi_c;
    
    ll(zz<=xyzor(3)) = ll_d;
    pp(zz<=xyzor(3)) = pp_d;
    rr(zz<=xyzor(3)) = rr_d;
    
    ll(zz>xyzor(3))  = ll_c;
    pp(zz>xyzor(3))  = pp_c;
    rr(zz>xyzor(3))  = rr_c;
end