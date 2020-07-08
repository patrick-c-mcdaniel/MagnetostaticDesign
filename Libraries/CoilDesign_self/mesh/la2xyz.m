function [ xout yout zout ] = la2xyz(lls, aas, r_in, r_out, f_cone, lst, lmid, lend )

    xout = zeros(size(lls));
    yout = zeros(size(lls));
    zout = zeros(size(lls));
    
    zout = zout + (lls>lmid).*lls + (lls<=lmid).*(2*lmid-lls);
    max(zout)
    lcon = lmid + (lend-lmid)*(f_cone);
    
    rrs = r_in * (lls<=lcon) + r_out * (lls>lmid) + (r_in + (r_out-r_in)*(lls-lcon)/(lmid-lcon)).*(lls>lcon).*(lls<=lmid);
    
    xout = cos(aas).*rrs;
    yout = sin(aas).*rrs;
 
end