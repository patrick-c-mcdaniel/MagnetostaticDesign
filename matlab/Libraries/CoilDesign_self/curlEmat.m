function [ curlEx, curlEy, curlEz ] = curlEmat( node,triangle,basis,r, w1  )

    [Cx,Cy,Cz] = Cn7_par(node,triangle,basis,r);
    
    curlEx = w1*Cx;
    curlEy = w1*Cy;
    curlEz = w1*Cz;

end