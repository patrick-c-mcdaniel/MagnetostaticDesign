function [ xout yout zout ] = la2xyz(lls, aas, r_in, r_out, f_cone, lst, lmid, lend )

    xout = zeros(size(lls));
    yout = zeros(size(lls));
    zout = zeros(size(lls));
    
    zout = (lls>lmid).*lls + ...
           (lls<=lmid).*(lls>(lmid-(r_out-r_in))).*lmid + ...
           (lls<=(lmid-(r_out-r_in))).*(lls+(r_out-r_in));
    max(zout)
    %lcon = lmid + (lend-lmid)*(f_cone);
    
    rrs =  (lls>lmid).*r_out + ...
           (lls<=lmid).*(lls>(lmid-(r_out-r_in))).*(lls-lmid+r_out) + ...
           (lls<=(lmid-(r_out-r_in))).*r_in;
    
    xout = cos(aas).*rrs;
    yout = sin(aas).*rrs;
 
end