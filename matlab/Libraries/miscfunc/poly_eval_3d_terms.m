function [ Vq ] = poly_eval_3d_terms( coefs3d, TERMS, xq, yq, zq )
%UNTITLED15 Summary of this function goes here
%   TERMS - (Nterms)x(3) array of polynomial coefficients
%               eg. [ 0 0 0; 1 0 0; 1 1 0; 1 0 1; 2 0 0 ] =>
%                       fits to (1), (x), (xy), (xz), (x2) terms


    Vq = zeros(size(xq));
    icoef = 1;
    
    Nter = size(TERMS,1);
    
    for iter = 1:Nter
        px = TERMS(iter, 1);
        py = TERMS(iter, 2);
        pz = TERMS(iter, 3);
        Vq = Vq + coefs3d(icoef)*(xq.^px).*(yq.^py).*(zq.^pz);
        icoef = icoef+1;

    end
    
end

