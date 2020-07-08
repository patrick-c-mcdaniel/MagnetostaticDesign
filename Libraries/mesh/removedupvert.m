function [ ver_red, fac_red, inds_unique ] = removedupvert( ver_all, fac_all, TOL )

    %% eliminate identical vertices


    map2min = zeros(size(ver_all,1),1);
    inds_unique = zeros(0, 1);
    ind_uni = 0;

    ver_red = zeros(0,3);
    fac_red = zeros(0,3);

    for ivv = 1:size(ver_all,1)

        indeq = find( sqrt( (ver_all(:,1)-ver_all(ivv,1)).^2 + (ver_all(:,2)-ver_all(ivv,2)).^2 + (ver_all(:,3)-ver_all(ivv,3)).^2 )<TOL );

        if ivv~=min(indeq)
            continue
        end
        ind_uni = ind_uni+1;

        ver_red = cat(1,ver_red, [ver_all(ivv,1), ver_all(ivv,2), ver_all(ivv,3)] );

        map2min(indeq) = ind_uni;
        inds_unique = cat(1, inds_unique, ivv);
    end


    for iff = 1:size(fac_all,1)


        fac_red = cat(1, fac_red, [map2min(fac_all(iff,1)), map2min(fac_all(iff,2)), map2min(fac_all(iff,3))] );
    end

end