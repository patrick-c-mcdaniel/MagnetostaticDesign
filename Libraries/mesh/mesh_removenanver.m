function [ ver_nonan, fac_nonan ] = mesh_removenanver( ver1, fac1 )

    Nver = size(ver1, 1);
    Nfac = size(fac1, 1);
    vernan = zeros(Nver,1);
    facnan = zeros(Nfac,1);
    
    for ivv = 1:Nver 
        if isnan(sum( ver1(ivv, :), 2))
            vernan(ivv) = 1;
        end
    end
    
    for iff = 1:Nfac 
        if sum( vernan( fac1(iff,:) ) )>0
            facnan(iff) = 1;
        end
    end

    fac_nonan = fac1( facnan==0, : );
    fac_nonan_orig = fac_nonan;
    ver_nonan = ver1( vernan==0, : );
    
    %% correct indices
    for ivv = 1:Nver
        if vernan(ivv)==0
            continue
        end
        fac_nonan( fac_nonan_orig > ivv ) = fac_nonan( fac_nonan_orig > ivv) - 1;
    end
end