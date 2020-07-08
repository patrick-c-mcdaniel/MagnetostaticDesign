% function [ cont_xyz ] = genisocont_sphcord( ver_mesh, fac_mesh, psi_mesh, xyz_ori, tz_bas, Niso, BSFN, PLOT, stdSM ) 
function [ cont_dec, len_cont ] = decimate_cont( varargin ) 
    %% function [cont_xyz] = genisocont_sphcord( ver_mesh, fac_mesh, psi_mesh, xyz_ori, tz_vec, Niso )
    %     Generate xyz isocontours from a stream function on a surface.
    %     This function parametrizes the surface by fitting it to spherical
    %     coordinates in the user-specified coordinate system. Operations
    %     are then performed on the "unfolded" {(theta, phi)} 2D surface
    %
    %   Input variables (required):
    %       ver_mesh        - (Nver)x(3) matrix of mesh vertices
    %       fac_mesh        - (Nfac)x(3) matrix of mesh triangles
    %       psi_mesh        - (Nver)x(1) vector of stream function values
    %       xyz_ori         - (1)x(3) vector specifying origin for
    %                           spherical coordinate system
    %       tz_bas         - (3)x(3) basis change matrix; rows are new
    %                           basis vectors
    %       Niso            - Number of isocontours to compute
    %   
    %   Output variables:
    %       cont_xyz        - cell array of (Ncont_i)x(3) matrices, with
    %                           each matrix containing a set of points
    %                           corresponding to 1 contour
    %
    %   Input variables (optional);
    %       BSFN            - Biot-Savart filename: include to write
    %                           biot-savart file for coil. Input [] will not write file
    %                           (Default: [])
    %       PLOT            - Plot isocontours in new figure (0 will not
    %                           plot; Default: 1)

    %% parse input
    p_in = inputParser;
    addRequired( p_in, 'contour');
    addRequired( p_in, 'DistTol');
    addRequired( p_in, 'AngTol');
%     addRequired( p_in, 'xyz_ori');
%     addRequired( p_in, 'tz_bas');
%     addRequired( p_in, 'Niso');
%     
%     addParameter( p_in, 'BSFN', []);
%     addParameter( p_in, 'PLOT', 1);
%     addParameter( p_in, 'stdSM', 2);
    
    parse(p_in,varargin{:});
    
    
    cont_xyz = p_in.Results.contour;
    DistTol = p_in.Results.DistTol;
    AngTol = p_in.Results.AngTol;
%     xyz_ori = p_in.Results.xyz_ori;
%     tz_bas = p_in.Results.tz_bas;
%     Niso = p_in.Results.Niso;
%     BSFN = p_in.Results.BSFN;
%     PLOT = p_in.Results.PLOT;
%     stdSM = p_in.Results.stdSM;
%     
    %% 
    Npts = size(cont_xyz,2);
    cont_dec = cont_xyz(:,1);
    cur_xyz = cont_xyz(:,1);
    len_cont = 0;
    
    for ipp = 1:(Npts-1)
        v12 = cont_xyz(:,ipp+1) - cur_xyz; 
        d12 = norm(v12);
        mask_d = d12>=DistTol;
        if size(cont_dec,2) < 2
            mask_a = 0;
        else
            v01 = cur_xyz - prev_xyz;
            d01 = norm(v01);
            n12 = v12/d12;
            n01 = v01/d01;
            c1  = n01' * n12;
            a1  = acos(c1);
            mask_a = a1>=AngTol;
        end
        if mask_d || mask_a
            cont_dec = cat(2, cont_dec, cont_xyz(:,ipp+1));
            
            prev_xyz = cur_xyz;
            cur_xyz = cont_xyz(:,ipp+1);
            
            len_cont = len_cont + norm(cur_xyz-prev_xyz);
            continue
        else
            continue
        end
        
    end
    
    figure(301); hold off
     plot3( cont_xyz(1,:), cont_xyz(2,:), cont_xyz(3,:) ); grid on; grid minor; axis equal;
    hold on;
     plot3( cont_dec(1,:), cont_dec(2,:), cont_dec(3,:) ); grid on; grid minor; axis equal;
    
end
