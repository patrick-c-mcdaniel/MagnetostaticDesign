function [ ver_mer, fac_mer ] = mergemeshes( ver1, fac1, ver2, fac2, TOL )

    ver_12 = cat(1, ver1, ver2);
    fac_12 = cat(1, fac1, fac2+size(ver1,1));
    [ver_red, fac_red] = removedupvert( ver_12, fac_12, TOL );
    [ver_mer, fac_mer] = removenullfac( ver_red, fac_red );
end