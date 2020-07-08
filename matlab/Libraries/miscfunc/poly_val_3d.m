function [ Vq ] = poly_val_3d( coefs3d, order, xq, yq, zq )

    Vq = zeros(size(xq));
    icoef = 1;
    
    for iordx = 0:order
        for iordy = 0:(order-iordx)
            for iordz = 0:(order-iordx-iordy)
                
                Vq = Vq + coefs3d(icoef)*(xq.^iordx).*(yq.^iordy).*(zq.^iordz);
                icoef = icoef+1;
            end
        end
    end
    
end

