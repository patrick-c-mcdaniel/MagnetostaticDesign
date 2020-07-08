function [ Gx Gy Gz ] = SH_grad_r( xx, yy, zz, l, m, del )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    
    rr = sqrt( xx.^2 + yy.^2 + zz.^2 );
    tt = acos( zz./rr );
    pp = atan2( yy, xx);
    
    [ sh ] = harmonicY_real( l, m, tt, pp,'type','complex','norm',true);
    sh_rr = sh .* (1./(rr.^(l+1)));
    
    %% x- derivative
    
    rrxp = sqrt( (xx+del/2).^2 + yy.^2 + zz.^2 );
    ttxp = acos( zz./rrxp );
    ppxp = atan2( yy, (xx+del/2));
    [ sh ] = harmonicY_real( l, m, ttxp, ppxp,'type','complex','norm',true);
    sh_xp = sh .* (1./(rrxp.^(l+1)));
    
    rrxm = sqrt( (xx-del/2).^2 + yy.^2 + zz.^2 );
    ttxm = acos( zz./rrxm );
    ppxm = atan2( yy, (xx-del/2));
    [ sh ] = harmonicY_real( l, m, ttxm, ppxm,'type','complex','norm',true);
    sh_xm = sh .* (1./(rrxm.^(l+1)));
    
    Gx = (sh_xp - sh_xm) / del;
    
    %% y-derivative
    
    rryp = sqrt( xx.^2 + (yy+del/2).^2 + zz.^2 );
    ttyp = acos( zz./rryp );
    ppyp = atan2( (yy+del/2), xx);
    [ sh ] = harmonicY_real( l, m, ttyp, ppyp,'type','complex','norm',true);
    sh_yp = sh .* (1./(rryp.^(l+1)));
    
    rrym = sqrt( xx.^2 + (yy-del/2).^2 + zz.^2 );
    ttym = acos( zz./rrym );
    ppym = atan2( (yy-del/2), xx);
    [ sh ] = harmonicY_real( l, m, ttym, ppym,'type','complex','norm',true);
    sh_ym = sh .* (1./(rrym.^(l+1)));
    
    Gy = (sh_yp - sh_ym) / del;
    
    %% z-derivative
    
    rrzp = sqrt( xx.^2 + yy.^2 + (zz+del/2).^2 );
    ttzp = acos( (zz+del/2)./rrzp );
    ppzp = atan2( yy, xx);
    [ sh ] = harmonicY_real( l, m, ttzp, ppzp,'type','complex','norm',true);
    sh_zp = sh .* (1./(rrzp.^(l+1)));
    
    rrzm = sqrt( xx.^2 + yy.^2 + (zz-del/2).^2 );
    ttzm = acos( (zz-del/2)./rrzm );
    ppzm = atan2( yy, xx);
    [ sh ] = harmonicY_real( l, m, ttzm, ppzm,'type','complex','norm',true);
    sh_zm = sh .* (1./(rrzm.^(l+1)));
    
    Gz = (sh_zp - sh_zm) / del;
end

