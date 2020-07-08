function [ coefs3d Vq ] = poly_fit_3d_terms( xx, yy, zz, VV, TERMS, xq, yq, zq )
%UNTITLED15 Summary of this function goes here
%   TERMS - (Nterms)x(3) array of polynomial coefficients
%               eg. [ 0 0 0; 1 0 0; 1 1 0; 1 0 1; 2 0 0 ] =>
%                       fits to (1), (x), (xy), (xz), (x2) terms

%     xx = 30*xx;
%     yy = 30*yy;
%     zz = 30*zz;
%     
%     
%     xq = 30*xq;
%     yq = 30*yq;
%     zq = 30*zq;
    Nter = size(TERMS,1);
    
    A = ones(numel(xx),0);
    for iter = 1:Nter
        px = TERMS(iter, 1);
        py = TERMS(iter, 2);
        pz = TERMS(iter, 3);
        A = [A (xx.^px).*(yy.^py).*(zz.^pz)];
    end
%     size(A)
%     size(dat1d)
    coefs3d = A\VV;
    
    %% evaluate at query points
    if nargin<6
        Vq = NaN;
        return
    end
    
    %%% v1 - memory intensive
%     Aq = ones(numel(xq),0);
% 
%     for iordx = 0:order
%         for iordy = 0:(order-iordx)
%             for iordz = 0:(order-iordx-iordy)
%                 
%                 Aq = [Aq (xq.^iordx).*(yq.^iordy).*(zq.^iordz)];
%             end
%         end
%     end
%     
%     Vq = Aq * coefs3d;
    
        %%% v2 - less memory intensive
    Vq = zeros(size(xq));
    icoef = 1;
    
    for iter = 1:Nter
        px = TERMS(iter, 1);
        py = TERMS(iter, 2);
        pz = TERMS(iter, 3);
        Vq = Vq + coefs3d(icoef)*(xq.^px).*(yq.^py).*(zq.^pz);
        icoef = icoef+1;

    end
    
end

