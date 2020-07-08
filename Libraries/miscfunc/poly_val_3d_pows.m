function [ Vq ] = poly_val_3d_pows( coefs3d, pows, xq, yq, zq )

    Vq = zeros(size(xq));
    
    
    for icoef = 1:numel(coefs3d)
        
        px = pow(icoef,1);
        py = pow(icoef,2);
        pz = pow(icoef,3);
        
        Vq = Vq + coefs3d(icoef)*(xq.^px).*(yq.^py).*(zq.^pz);

    end
    
end

