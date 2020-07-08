function [ ver_red fac_red ] = removenullfac( ver_all, fac_all )

    ver_red = ver_all;
    fac_red = zeros(0, 3);
    for iff = 1:size(fac_all,1)

        fac_tmp = fac_all(iff,:);
        if numel(unique(fac_tmp))<3
            continue;
        end
        fac_red = cat(1, fac_red, fac_tmp);
        
    end

end