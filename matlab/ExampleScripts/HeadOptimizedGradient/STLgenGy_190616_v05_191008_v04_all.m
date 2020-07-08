% 
% %% interp. vertices (at native mesh resolution)
% 
% [ ver_up, fac_up, psi_up ] = meshScalarUpsample( coil_red.listNode, coil_red.listTriangle, coil_red.s, 3 );
% 
% %% generate isocontours
% 
% 
% displayStreamFunction(fac_up,psi_up,ver_up);
% 
% ver_mesh = ver_up;
% fac_mesh = fac_up;
% psi_mesh = psi_up;
% xyz_ori = [ 0.15 0 0 ];
% tz_bas = [ 0 0 -1 ; 0 -1 0 ; -1 0 0 ];
% Niso = 30;
% 
% BSFN = 'BSfiles/coilGy_190614_v04p0.biot';
% 
% cont_Gy_all = genisocont_sphcord( ver_mesh, fac_mesh, psi_mesh, xyz_ori, tz_bas, Niso,'BSFN',[] );

%% copied from genMesh_Gy_190614_v04p0

% Parameters for the mesh, coil, and wire grooves

Rcform = 0.178889;
Rsform = 0.191299;
WGRV = 1.074e-3;
TGRV = WGRV*1.05;
TBOT = 1.5e-3;

Dcyl = 2*Rcform + TBOT;
Dsph = 2*Rsform + TBOT;
Dfpl = Dcyl - 2*0.04;

Rcyl = Dcyl/2;
Rsph = Dsph/2;
Rfpl = Dfpl/2;

Ncir = 2^6;

Ccyl = pi*Dcyl;
Csph = pi*Dsph;

delR_e = Ccyl/Ncir/2;

Xmax = 0.177; % + 0.04;
Xbnd = 0; %( Rsph^2 - Rcyl^2 )^(1/2);

Lcyl = Xmax-Xbnd;

Nlcyl = ceil(Lcyl/delR_e);

Nfpl = 8;

delL = Lcyl/Nlcyl;

RATIO_RC = 0.7;

%% look at contours; generate normal vectors

load('./optresults/cont_Gy_all_Example_200626.mat');

figure(201); hold off;
figure(202); hold off;
Npt_all = 0;
% cont_Gy_all = cont_Gx_all;

cont_Gy_all_dec = cell(size(cont_Gy_all));  % cell array containing all decimated contours
rn_all_dec = cell(size(cont_Gy_all));       % cell array containing surface normal vectors for all decimated contour points
                                            %  This is computed for the mesh geometry; knowing this helps determine the direction of the
                                            %  curve, the direction along which to "thicken" the curve, and the direction into the
                                            %  volumetric coil former

len_all = 0;

TOL_end = 1e-6;
TolXmax = 2e-3;

for iccc = 1:1:numel(cont_Gy_all)
    cont_Gy_1_undec = cont_Gy_all{iccc};
    
    DistTol = 7e-3;
    AngTol  = pi/180*360;
    
    [ cont_Gy_1, len_cont ] = decimate_cont( cont_Gy_1_undec, DistTol, AngTol );  % decimates contour based on length/angle tolerances
    if norm( cont_Gy_1(:,1)-cont_Gy_1(:,end) )<TOL_end
        cont_Gy_1 = cont_Gy_1(:,1:end-1);
    end
    
    len_all = len_all + len_cont;
    
figure(201);
    plot3( cont_Gy_1(1,:), cont_Gy_1(2,:), cont_Gy_1(3,:) ); grid on; grid minor; axis equal;
    hold on;
    figure(202);
    plot3( cont_Gy_1(1,:), cont_Gy_1(2,:), cont_Gy_1(3,:) ); grid on; grid minor; axis equal;
    hold on;
    %%% compute normal vectors

    x_1 = cont_Gy_1(1,:);
    y_1 = cont_Gy_1(2,:);
    z_1 = cont_Gy_1(3,:);

    r_1  = zeros(size(cont_Gy_1));
    rn_1 = zeros(size(cont_Gy_1));


    msph = x_1<Xbnd;
    mcyl = x_1>=Xbnd & x_1<(Xmax-TolXmax);
    max(x_1(:))
    mfil = x_1>=(Xmax-TolXmax);
    mann = x_1>=(Xmax-TOL_end);
    %%% cyl normals
    r_1(:,mcyl) = cat(1, zeros(size(y_1(mcyl))), y_1(mcyl), z_1(mcyl));
    r_1(:,msph) = cat(1, x_1(msph), y_1(msph), z_1(msph));
    r_xy = sqrt(y_1.^2 + z_1.^2);
    r_1(:,mfil) = cat(1, (x_1(mfil)-Xmax+TolXmax), y_1(mfil)./ r_xy(mfil) .* sqrt( TolXmax^2-(x_1(mfil)-Xmax+TolXmax).^2), z_1(mfil)./ r_xy(mfil) .* sqrt( TolXmax^2-(x_1(mfil)-Xmax+TolXmax).^2));
    r_1(:,mann) = cat(1, ones(size(x_1(mann))), zeros(size(x_1(mann))), zeros(size(x_1(mann))));
    
    mr_1 = sqrt( sum( r_1.^2, 1));
    
    rn_1 = r_1 ./ repmat(mr_1, [3 1]);
figure(201);
    quiver3( x_1, y_1, z_1, rn_1(1,:), rn_1(2,:), rn_1(3,:) );
    Npt = size(cont_Gy_1,2);
    Npt_all = Npt + Npt_all;
    
    cont_Gy_all_dec{iccc} = cont_Gy_1;
    rn_all_dec{iccc} = rn_1;
end

%% compute on-surface normal vectors 

% rn_all_dec : cell array of surface normals
% cont_Gy_all_dec : cell array of contours (points)

v_pos_all = cell( size( rn_all_dec ) );


for iccc = 1:1:numel(cont_Gy_all_dec)
    c_tmp = cont_Gy_all_dec{iccc};
    if isempty(c_tmp)
        continue;
    end
    rn_tmp = rn_all_dec{iccc};
    
    v_c_tmp = zeros(size(c_tmp));
    v_c_tmp(:,1:end-1) = c_tmp(:,2:end)-c_tmp(:,1:end-1);
    v_c_tmp(:,end) = c_tmp(:,1)-c_tmp(:,end);
    
    m_c_tmp = sqrt( sum(v_c_tmp.^2,1) );
    n_c_tmp = v_c_tmp ./ repmat(m_c_tmp,[3 1]);
    figure(401); hold off;
    plot3( c_tmp(1,:), c_tmp(2,:), c_tmp(3,:),'o' ); grid on; grid minor; axis equal;
    hold on;
    quiver3( c_tmp(1,:), c_tmp(2,:),c_tmp(3,:), n_c_tmp(1,:), n_c_tmp(2,:), n_c_tmp(3,:) ); 
    
    
    %%% create perpendicular vectors 
    v_pc_tmp = cross( rn_tmp, n_c_tmp, 1);
    m_pc_tmp = sqrt( sum( v_pc_tmp.^2, 1));
    n_pc_tmp = v_pc_tmp ./ repmat( m_pc_tmp,[3 1]);
    ca_tmp = zeros(1, size(c_tmp,2));
    ca_tmp(:,2:end) = sqrt( 1/2*(1+ diag(n_pc_tmp(:,2:end)'*n_pc_tmp(:,1:end-1))));
    ca_tmp(:,1) = sqrt( 1/2*(1+ n_pc_tmp(:,1)'*n_pc_tmp(:,end)));
    
    v12 = ( n_pc_tmp(:,1:end) + circshift( n_pc_tmp, 1, 2 ));
    m12 = sqrt( sum( v12.^2, 1 ));
    n12 = v12 ./ repmat(m12, [3 1]);
    
    v_pos = repmat(1./ca_tmp, [3 1]) .*  n12;
    
    figure(501); hold off;
    plot3( c_tmp(1,:), c_tmp(2,:), c_tmp(3,:),'o' ); grid on; grid minor; axis equal;
    
    hold on;
    quiver3( c_tmp(1,:), c_tmp(2,:),c_tmp(3,:), v_pos(1,:), v_pos(2,:), v_pos(3,:) ); 
    
    v_pos_all{iccc} = v_pos;
end

%% create inner/outer windings for on-surface contours
%     These will represent the inside and outside edges of the thickened contours 

c_in_all  = cell( size( rn_all_dec ));
c_out_all = cell( size( rn_all_dec ));

ver_all = cell( size( rn_all_dec ));
fac_all = cell( size( rn_all_dec ));
ver_all3d = cell( size( rn_all_dec ));
fac_all3d = cell( size( rn_all_dec ));


% WGRV = 1.074e-3;
TGRV_etch = 20e-3;

INDBACK = 1:83;
% INDBACK = [ 1:6 8 10 12:17 20:22 24:25 27 28 30 32:34 37:48 50:70 72:83 ];
% INDBACK = [ 7 9 11 18 19 23 26 29 31 35 36 43 49 71 ];

figure(701); hold off;

for iccc = 1:1:numel(cont_Gy_all_dec)
    c_tmp = cont_Gy_all_dec{iccc};
    if isempty(c_tmp)
        continue;
    end
    v_pos_tmp = v_pos_all{iccc};
    rn_tmp = rn_all_dec{iccc};
    
    c_in_tmp = c_tmp + WGRV/2*v_pos_tmp;
    c_out_tmp = c_tmp - WGRV/2*v_pos_tmp;
    
    figure(601); hold off;
    plot3( c_in_tmp(1,:), c_in_tmp(2,:), c_in_tmp(3,:) ); grid on; grid minor; axis equal;
    hold on;
    plot3( c_out_tmp(1,:), c_out_tmp(2,:), c_out_tmp(3,:) ); grid on; grid minor; axis equal;
   
    c_in_all{iccc} = c_in_tmp;
    c_out_all{iccc} = c_out_tmp;
    
    Ntmp = size(c_tmp,2);
    
    
    tri_o = cat(1, 1:Ntmp, (Ntmp+1):2*Ntmp, [2:Ntmp, 1]);
    tri_i = cat(1, (Ntmp+1):2*Ntmp, [(Ntmp+2):2*Ntmp, Ntmp+1], [2:Ntmp, 1]);
    
    %%% further triangles for 3D mesh
    
    
    tri_bc1 = cat(1, 1:Ntmp, (Ntmp+1):2*Ntmp, [(Ntmp+2):2*Ntmp, Ntmp+1]) + Ntmp;
    tri_bc2 = cat(1, 1:Ntmp, [(Ntmp+2):2*Ntmp, Ntmp+1], [2:Ntmp, 1]) + Ntmp;
    tri_cd1 = tri_o + 2*Ntmp;
    tri_cd2 = tri_i + 2*Ntmp;
    tri_da1 = tri_o;
    tri_da2 = tri_i;
    tri_da1([1 3],:) = tri_da1([1 3],:) +3*Ntmp;
    tri_da1([2],:) = tri_da1(2,:) - Ntmp;
    tri_da2(3,:) = tri_da2(3,:) + 3*Ntmp;
    tri_da2([1 2],:) = tri_da2([1 2],:) - Ntmp;
    
    ver_tmp = [ c_out_tmp, c_in_tmp ];
    fac_tmp = [ tri_o, tri_i ];
    ver_all{iccc} = ver_tmp;
    fac_all{iccc} = fac_tmp;
    
    %%% ver/fac for 3D mesh
    ver_tmp3d = [ c_out_tmp, c_in_tmp, c_in_tmp+TGRV_etch*rn_tmp, c_out_tmp+TGRV_etch*rn_tmp ];
    fac_tmp3d = [ tri_o, tri_i, tri_bc1, tri_bc2, tri_cd1, tri_cd2, tri_da1, tri_da2 ];
    ver_all3d{iccc} = ver_tmp3d;
    fac_all3d{iccc} = fac_tmp3d;
    
    figure(701); 
    patch('Faces',fac_tmp','Vertices',ver_tmp','FaceColor',[1, 0.8, 0.8],'EdgeColor',[0 0 0]);
    hold on;
    
    figure(801); 
    patch('Faces',fac_tmp3d','Vertices',ver_tmp3d','FaceColor',[1, 0.8, 0.8],'EdgeColor',[0 0 0]);
    hold on;
    
end

%% write test 2D stl
%      This is an STL of the 2D coil - the 1D contours have been turned into 2D meshes, corresponding to the projections of the coil wire
%      windings onto the coil former surface.

ver_all_mat = zeros(3, 0);
fac_all_mat = zeros(3, 0);
Ncum = 0;

for iccc = 1:1:numel(cont_Gy_all_dec)
    ver_tmp = ver_all{iccc};
    if isempty(ver_tmp)
        continue;
    end
    fac_tmp = fac_all{iccc};
    Ntmp = size(fac_tmp,2);
    
    ver_all_mat = [ver_all_mat, ver_tmp];
    fac_all_mat = [fac_all_mat, fac_tmp+Ncum];
    
    Ncum = Ncum + Ntmp;
    
end

FN_stl = './STL/STL2D/Gy2D_190916_v05_191008_v04.stl';
stlwrite( FN_stl, fac_all_mat', ver_all_mat' );

%% write 3D STLs for each individual winding
%       These can then be combined in a mesh editor/CAD program (such as inventor)
%       This avoids the problem of generating a single 3D STL file with intersecting windings (and thus, intersecting meshes). The
%       intersections can be dealth with when merging the mesh-derived volumes in Inventor, for example.

ver_all3d_mat = zeros(3, 0);
fac_all3d_mat = zeros(3, 0);
Ncum = 0;

for iccc = 1:1:numel(cont_Gy_all_dec)

    ver_tmp = ver_all3d{iccc};
    fac_tmp = fac_all3d{iccc};
    Ntmp = size(ver_tmp,2);
        if isempty(ver_tmp)
        continue;
    end
    ver_all3d_mat = [ver_all3d_mat, ver_tmp];
    fac_all3d_mat = [fac_all3d_mat, fac_tmp+Ncum];
    
    Ncum = Ncum + Ntmp;
    
    FN_stl = ['./STL/Gy_sep/Gy_191008_v05_3d_sep_v04' num2str(iccc,'%02u') '.stl'];
    stlwrite( FN_stl, fac_tmp', ver_tmp' );

end

    %% Generate 3D stl with all wire windings
    
    FN_stl = ['./STL/STL3D/Gy_191008_v05p4_3d.stl'];
    stlwrite( FN_stl, fac_all3d_mat', ver_all3d_mat' )
%% combine RO/shim coil with magnet model in single BS-file

% BSFN_mag = '/home/patmcd/Documents/HelmetMagnet/Magnet_v2/190327_proc/BSfiles/Magnet_190327_v01.biot';
% 
% BSFN_mag_ShimRO = 'BSfiles/Mag_ShimROy_190402_v04_test.biot';
% 
% BSfile_combine( BSFN_mag, BSFN, BSFN_mag_ShimRO );
