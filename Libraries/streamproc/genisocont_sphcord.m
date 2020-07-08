% function [ cont_xyz ] = genisocont_sphcord( ver_mesh, fac_mesh, psi_mesh, xyz_ori, tz_bas, Niso, BSFN, PLOT, stdSM ) 
function [ cont_xyz ] = genisocont_sphcord( varargin ) 
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
    addRequired( p_in, 'ver_mesh');
    addRequired( p_in, 'fac_mesh');
    addRequired( p_in, 'psi_mesh');
    addRequired( p_in, 'xyz_ori');
    addRequired( p_in, 'tz_bas');
    addRequired( p_in, 'Niso');
    
    addParameter( p_in, 'BSFN', []);
    addParameter( p_in, 'PLOT', 1);
    addParameter( p_in, 'stdSM', 2);
    
    parse(p_in,varargin{:});
    
    
    ver_mesh = p_in.Results.ver_mesh;
    fac_mesh = p_in.Results.fac_mesh;
    psi_mesh = p_in.Results.psi_mesh;
    xyz_ori = p_in.Results.xyz_ori;
    tz_bas = p_in.Results.tz_bas;
    Niso = p_in.Results.Niso;
    BSFN = p_in.Results.BSFN;
    PLOT = p_in.Results.PLOT;
    stdSM = p_in.Results.stdSM;
    
    load('cm_rb.mat');
    
    tz_bas_inv = inv(tz_bas);
    %% generate 2D surface
    
    x3d = ver_mesh(:,1);
    y3d = ver_mesh(:,2);
    z3d = ver_mesh(:,3);
    
    x3d_or = x3d - xyz_ori(1);
    y3d_or = y3d - xyz_ori(2);
    z3d_or = z3d - xyz_ori(3);
    
    xyz3_or = cat(1, x3d_or(:)', y3d_or(:)', z3d_or(:)' );
    xyz3_bas = tz_bas * xyz3_or;
    
    x3_bas = xyz3_bas(1,:);
    y3_bas = xyz3_bas(2,:);
    z3_bas = xyz3_bas(3,:);
    
    rr3 = sqrt( x3_bas.^2 + y3_bas.^2 + z3_bas.^2 );
    tt3 = acos( z3_bas./rr3);
    pp3 = atan2( y3_bas, x3_bas );
    
    %% generate interpolation surface
    
    thtint = linspace(min(tt3(:)), max(tt3(:)), 2048);
    phiint = linspace(-pi, pi, 2048);

    [ thtint2d phiint2d ] = ndgrid( thtint, phiint );
%     
%     Int_psi = scatteredInterpolant( tt3(:), pp3(:), psi_mesh(:),'natural');
%     Int_rr  = scatteredInterpolant( tt3(:), pp3(:), rr3(:),'natural');


    Int_psi = scatteredInterpolant( tt3(:).*cos(pp3(:)), tt3(:).*sin(pp3(:)), psi_mesh(:),'linear');
    Int_rr  = scatteredInterpolant( tt3(:).*cos(pp3(:)), tt3(:).*sin(pp3(:)), rr3(:),'linear');

    rrint   =  Int_rr(  thtint2d(:).*cos(phiint2d(:)), thtint2d(:).*sin(phiint2d(:)) );
    psiint  =  Int_psi( thtint2d(:).*cos(phiint2d(:)), thtint2d(:).*sin(phiint2d(:)) );
    
    %%% "x-y cartesian" interpolate
    
    xtp = linspace( min(thtint2d(:).*cos(phiint2d(:))), max(thtint2d(:).*cos(phiint2d(:))), 2048 );
    ytp = linspace( min(thtint2d(:).*sin(phiint2d(:))), max(thtint2d(:).*sin(phiint2d(:))), 2048 );
    
    [ x2tp y2tp ] = ndgrid( xtp, ytp );
    
    rrxy   =  Int_rr(  x2tp, y2tp );
    psixy  =  Int_psi( x2tp, y2tp );
    
    rtp_max = max( thtint2d(:) );
    
    psixy( (x2tp.^2 + y2tp.^2) > rtp_max^2) = 0;

    zzint   = rrint(:) .* cos( thtint2d(:) );
    yyint   = rrint(:) .* sin( thtint2d(:) ) .* sin( phiint2d(:) );
    xxint   = rrint(:) .* sin( thtint2d(:) ) .* cos( phiint2d(:) );

    figure(9); hold off;
    
    cdat_all = dat2colors( psiint(:), cm_rb);
    
    scatter( thtint2d(:).*cos(phiint2d(:)), thtint2d(:).*sin(phiint2d(:)), 10, cdat_all ); %colormap(cm_rb); 
    axis equal
%     surf( thtint2d, phiint2d, psiint,'Edgecolor','none' ); colormap jet; view(2)
    

    %% smooth interp stream function
    psixy_smth = imgaussfilt(psixy, stdSM,'Padding',0);
    psixy_smth( (x2tp.^2 + y2tp.^2) > rtp_max^2) = 0;
    
%     rrxy_smth = 
    %% generate isocontours
    figure(10);
    imagesc( psixy_smth ); axis image; colormap(cm_rb);
    
    figure(11); 
    [ cont_all ] = contour( x2tp, y2tp, psixy_smth, Niso ); colormap(cm_rb); axis equal;
    
    cont_xy_all = cell(1,0);
    
    Ncont = 0;
    icc_tmp = 1;
    
    ind_ca = 1;
    %%% separate individual contours
    while ind_ca<size(cont_all,2)
        cont_xy_tmp = cont_all(:, (ind_ca+1):(ind_ca+cont_all(2, ind_ca)));
        cont_xy_all{icc_tmp} = cont_xy_tmp;
        ind_ca = ind_ca + cont_all(2, ind_ca)+1;
        Ncont = Ncont + 1;
        icc_tmp = icc_tmp+1;
    end
    
    cont_xyz_all = cell(1, Ncont);
    cont_xyz_ori = cell(1, Ncont);
    
    %%% convert contours to (theta, phi) coordinates; then (x,y,z)
    %%% coordinates
    
    for icc = 1:Ncont
        cont_tmp = cont_xy_all{icc};
        
        rr_tmp   =  Int_rr(  cont_tmp(1,:), cont_tmp(2,:) );
        
        tt_tmp = sqrt( cont_tmp(1,:).^2 + cont_tmp(2,:).^2 );
        pp_tmp = atan2(  cont_tmp(2,:), cont_tmp(1,:) );
        
        cont_tp_all = [ tt_tmp(:), pp_tmp(:)]';
        
        xx_tmp = rr_tmp .* sin(tt_tmp) .* cos(pp_tmp);
        yy_tmp = rr_tmp .* sin(tt_tmp) .* sin(pp_tmp);
        zz_tmp = rr_tmp .* cos(tt_tmp);
        
        cont_xyz_all{icc} = [ xx_tmp(:) yy_tmp(:) zz_tmp(:) ]';
        
        cont_xyz_ori_tmp = tz_bas_inv * cont_xyz_all{icc};
        cont_xyz_ori_tmp(1,:) = cont_xyz_ori_tmp(1,:) + xyz_ori(1);
        cont_xyz_ori_tmp(2,:) = +1*( cont_xyz_ori_tmp(2,:) + xyz_ori(2) );
        cont_xyz_ori_tmp(3,:) = +1*( cont_xyz_ori_tmp(3,:) + xyz_ori(3) );
        cont_xyz_ori{icc} = cont_xyz_ori_tmp;
    end
    
    %%% plot isocontours
    
    figure(12); hold off;
    for icc = 1:Ncont
        cont_tmp = cont_xyz_ori{icc};
        plot3( cont_tmp(1,:), cont_tmp(2,:), cont_tmp(3,:) ); grid on; grid minor; axis equal;
        hold on;
    end
    
    %% map isocontours back to original surface
    
    
    %%% output variable
    
    cont_xyz = cont_xyz_ori;
    

    %% write BS file
    
    if isempty(BSFN)
        return
    end
    
    writeBScont_xyz(cont_xyz, BSFN);
    
end
