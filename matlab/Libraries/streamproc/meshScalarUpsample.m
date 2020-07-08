function [ ver_up, fac_up, psi_up ] = meshScalarUpsample( ver, fac, psi, Nup )

    Nfac = size(fac,1);
    Nver = size(ver,1);
    
    Nver_up = Nver * ((Nup+1)*Nup/2);
    Nfac_up = Nfac * Nup^2;
    
    ver_up = zeros(Nver_up, 3);
    fac_up = zeros(Nfac_up, 3);
    psi_up = zeros(Nver_up, 1);
    
    %% matrices for interp/new triangles
    
    Avup = zeros( (Nup+1)*Nup/2, 3 );
    
    irr = 1;
    for in1 = 1:(Nup+1)
        for in2 = 1:in1
               Avup( irr,: ) = [ (Nup+1)-in1, in1-in2, in2-1 ];
               irr = irr+1;
        end
    end
    Avup = 1/Nup*Avup;
    
    Fmat = zeros( Nup^2, 3 );
    irr = 1;
    
    for in1 = 1:(Nup)
        tri  = (in1)*(in1-1)/2;
        trip = (in1)*(in1+1)/2;
        tripp = (in1+1)*(in1+2)/2;
        
        for in2 = 1:(in1-1)
            Fmat(irr,:) = [ 1+tri+(in2-1), 1+trip+(in2-1), 1+trip+(in2)  ]; 
            irr = irr+1;
            Fmat(irr,:) = [ 1+tri+(in2-1), 1+trip+(in2),   1+tri+(in2)  ]; 
            irr = irr+1;
        end
        Fmat(irr,:) = [ trip, tripp-1, tripp ];
        irr = irr+1;
    end
    
    
    %% get net vertices/scalars/faces
    ivv = 1;
    ifc = 1;
    for iff = 1:Nfac
        
        fac_up(ifc:(Nup^2+ifc-1),:) = ivv-1 + Fmat;
        ifc = ifc+Nup^2;
        
        v1 = ver(fac(iff,1),:);
        v2 = ver(fac(iff,2),:);
        v3 = ver(fac(iff,3),:);
        v123 = cat(1,v1,v2,v3);
        vnew = Avup * v123;
        ver_up(ivv:(ivv+(Nup+2)*(Nup+1)/2-1),:) = vnew;
       
        
        p1 = psi(fac(iff,1));
        p2 = psi(fac(iff,2));
        p3 = psi(fac(iff,3));
        p123 = cat(1,p1,p2,p3);
        pnew = Avup * p123;
        psi_up(ivv:(ivv+(Nup+2)*(Nup+1)/2-1)) = pnew;
        ivv = ivv+(Nup+2)*(Nup+1)/2;
        
        
    end
    
    [ver_dv, fac_dv, inds_unique] = removedupvert( ver_up, fac_up, 1e-6 );
    [ver_nf, fac_nf] = removenullfac( ver_dv, fac_dv);
    ver_up = ver_nf;
    fac_up = fac_nf;
    
    psi_dv = psi_up(inds_unique);
    psi_up = psi_dv;
    
    figure(101); patch('Faces',fac_up,'Vertices',ver_up,'FaceColor',[1 0.7 0.7]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal

end