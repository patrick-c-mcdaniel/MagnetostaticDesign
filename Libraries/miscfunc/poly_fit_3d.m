function [ coefs3d, Vq, powmat ] = poly_fit_3d( xx, yy, zz, VV, order, xq, yq, zq )
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here

%     xx = 30*xx;
%     yy = 30*yy;
%     zz = 30*zz;
%     
%     
%     xq = 30*xq;
%     yq = 30*yq;
%     zq = 30*zq;
    
    A = ones(numel(xx),0);
    powmat = zeros( 0, 3 );
    
    for iordx = 0:order
        for iordy = 0:(order-iordx)
            for iordz = 0:(order-iordx-iordy)
                
                A = [A (xx(:).^iordx).*(yy(:).^iordy).*(zz(:).^iordz)];
                powmat = cat(1, powmat, [iordx, iordy, iordz]);
            end
        end
    end
%     size(A)
%     size(dat1d)
    coefs3d = A\VV(:);
    
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
    
    for iordx = 0:order
        for iordy = 0:(order-iordx)
            for iordz = 0:(order-iordx-iordy)
                
                Vq = Vq + coefs3d(icoef)*(xq.^iordx).*(yq.^iordy).*(zq.^iordz);
                icoef = icoef+1;
            end
        end
    end
    
end

