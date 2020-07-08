function [ Mxyz, AMred, A_comp ] = Mred2Mxyz_A_comp_4quad( Mred, inds_pp, inds_pn, inds_nn, inds_np, Lhigh )
% converts full Mxyz vector to reduced Mred vector. Mred is a compressed,
% symmetrized description of Mxyz; use of Mred imposes symmetry on the
% magnet design and reduces the size of the field computation/optimization
% problem
%
% Inputs:
%       Mred    - Reduced magnet vector
%       Imid    - Indices of magnets (into Mxyz) (pptx convention) located along the
%                   XZ-plane (symmetry plane). Keeps x- and z- components
%                   for these blocks; y-components forced to be ==0;
%       Iypos   - Indices of blocks with positive y-coordinate (must be
%                   sorted in ascending order). 
%       Iyneg   - Indices of blocks with negative y-coordinate (must
%                   pairwise match Iypos)
%       Lhigh   - set of spherical harmonic orders >1 (for generating Aexp
%                   matrix size; must be in increasing order); default=[]
%
% Outputs:
%       Mxyz    - Full magnet vector (eg for input into field computation)
%       Aexp    - Expansion matrix (ie: Mxyz = Aexp * Mred)
%       Acomp   - Expansion matrix for field computation (ie: D135red =
%                   D135*Acomp);
    
    if nargin<5
        Lhigh = [];
    end
    
    Npp = numel(inds_pp);
    Npn = numel(inds_pn);
    Nnn = numel(inds_nn);
    Nnp = numel(inds_np);
    
%     if Npos~=Nneg
%         error('Iypos and Iyneg must have same number of elements');
%     end
    Mxyz = zeros(3*(Npp+Npn+Nnn+Nnp),1);
    
    Mxyz(3*(inds_pp-1)+1) =  Mred(1:3:(3*Npp));
    Mxyz(3*(inds_pp-1)+2) =  Mred(2:3:(3*Npp));
    Mxyz(3*(inds_pp-1)+3) =  Mred(3:3:(3*Npp));
    
    Mxyz(3*(inds_pn-1)+1) = -Mred(1:3:(3*Npp));
    Mxyz(3*(inds_pn-1)+2) = -Mred(2:3:(3*Npp));
    Mxyz(3*(inds_pn-1)+3) =  Mred(3:3:(3*Npp));
    
    Mxyz(3*(inds_nn-1)+1) = -Mred(1:3:(3*Npp));
    Mxyz(3*(inds_nn-1)+2) =  Mred(2:3:(3*Npp));
    Mxyz(3*(inds_nn-1)+3) =  Mred(3:3:(3*Npp));
    
    Mxyz(3*(inds_np-1)+1) =  Mred(1:3:(3*Npp));
    Mxyz(3*(inds_np-1)+2) = -Mred(2:3:(3*Npp));
    Mxyz(3*(inds_np-1)+3) =  Mred(3:3:(3*Npp));
    
    %% create Aexp matrix (use to compress D-matrix)
    AMred = zeros( numel(Mxyz), numel(Mred) );
    
    Nbred = numel(Imid) + numel(Iypos);
    
    A_comp = zeros( numel(Mxyz), Nbred*3);
    
    for iim = 1:numel(Imid)
       AMred( 3*(Imid(iim)-1)+1, 1+2*(iim-1) ) = 1;
       AMred( 3*(Imid(iim)-1)+3, 2+2*(iim-1) ) = 1;
       
       A_comp( 3*(Imid(iim)-1)+1, 1+3*(iim-1) ) = 1;
       A_comp( 3*(Imid(iim)-1)+3, 3+3*(iim-1) ) = 1;
       
    end
    for iip = 1:numel(Iypos) 
       AMred( 3*(Iypos(iip)-1)+1, (2*Nmid+1)+3*(iip-1) ) = 1;
       AMred( 3*(Iypos(iip)-1)+2, (2*Nmid+2)+3*(iip-1) ) = 1;
       AMred( 3*(Iypos(iip)-1)+3, (2*Nmid+3)+3*(iip-1) ) = 1;
       
       A_comp( 3*(Iypos(iip)-1)+1, (3*Nmid+1)+3*(iip-1) ) = 1;
       A_comp( 3*(Iypos(iip)-1)+2, (3*Nmid+2)+3*(iip-1) ) = 1;
       A_comp( 3*(Iypos(iip)-1)+3, (3*Nmid+3)+3*(iip-1) ) = 1;
    end
    for iin = 1:numel(Iyneg)
       AMred( 3*(Iyneg(iin)-1)+1, (2*Nmid+1)+3*(iin-1) ) = 1;
       AMred( 3*(Iyneg(iin)-1)+2, (2*Nmid+2)+3*(iin-1) ) = -1;
       AMred( 3*(Iyneg(iin)-1)+3, (2*Nmid+3)+3*(iin-1) ) = 1;
       
       A_comp( 3*(Iyneg(iin)-1)+1, (3*Nmid+1)+3*(iin-1) ) = 1;
       A_comp( 3*(Iyneg(iin)-1)+2, (3*Nmid+2)+3*(iin-1) ) = -1;
       A_comp( 3*(Iyneg(iin)-1)+3, (3*Nmid+3)+3*(iin-1) ) = 1;
    end
    
    if numel(Lhigh)==0
%         AMred = Aexp1;
        return
    end
    
    %% go through higher-order SH terms; make bigger A-matrix; 
    %       Note: L=1 is dealt with separately as the order for L=1 terms
    %       (dipole terms) is different from higher-order terms. L=1 terms
    %       are ordered in "x-y-z" order, which is not the "m=-1,0,+1"
    %       order. Higher-order terms are ordered in increasing order of m
    
    Nbk = numel(Mxyz)/3;
    
    for ill = 1:numel(Lhigh)
       
        ll = Lhigh(ill);
        Nm = 2*ll+1;
        Mset = -ll:ll;
        
        A_ll = zeros( Nbk*Nm, Nbred*Nm );
        
        for imm = 1:Nm
            mm = Mset(imm);
            
            for iim = 1:numel(Imid)
               
                A_ll( Nm*(Imid(iim)-1)+imm, imm+Nm*(iim-1) ) = 0*(mm<0) + 1*(mm>=0);
                
            end
            
            for iip = 1:numel(Iypos) 
                A_ll( Nm*(Iypos(iip)-1)+imm, (Nm*Nmid+imm)+Nm*(iip-1) ) = 1*(mm<0) + 1*(mm>=0);
                
            end
            for iin = 1:numel(Iyneg) 
                A_ll( Nm*(Iyneg(iin)-1)+imm, (Nm*Nmid+imm)+Nm*(iin-1) ) = -1*(mm<0) + 1*(mm>=0);
                
            end
        end
                
        A_comp = [ A_comp,                               zeros(size(A_comp,1), size(A_ll,2)) ; ...
                   zeros(size(A_ll,1),size(A_comp,2)),   A_ll                                ];
        
        
    end
    
%     
%     Aexp( 3*(Imid-1)+1, 1:2:(2*Nmid) ) = 1;
%     Aexp( 3*(Imid-1)+3, 2:2:(2*Nmid) ) = 1;
    
%     Aexp( 3*(Iypos-1)+1, (2*Nmid+1):3:end ) = 1;
%     Aexp( 3*(Iyneg-1)+1, (2*Nmid+1):3:end ) = 1;
%     Aexp( 3*(Iypos-1)+2, (2*Nmid+2):3:end ) = 1;
%     Aexp( 3*(Iyneg-1)+2, (2*Nmid+2):3:end ) = -1;
%     Aexp( 3*(Iypos-1)+3, (2*Nmid+3):3:end ) = 1;
%     Aexp( 3*(Iyneg-1)+3, (2*Nmid+3):3:end ) = 1;  
    
%     Mxyz - Aexp*Mred

end

