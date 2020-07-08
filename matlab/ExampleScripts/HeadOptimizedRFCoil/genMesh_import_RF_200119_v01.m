
Dcyl = 0.205;
Dsph = 0.205;

Rcyl = Dcyl/2;
Rsph = Dsph/2;

Ncir = 2^6;

Ccyl = pi*Dcyl;
Csph = pi*Dsph;

delR_e = Ccyl/Ncir/2;

Lcyl = 0.168;

Xmax = 0.168;

Nlcyl = ceil(Lcyl/delR_e);

delL = Lcyl/Nlcyl;

RATIO_RC = 0.7;



%% create cylinder mesh

[ ver_cyl, fac_cyl ] = gencylmeshz( Xmax, Xmax-Lcyl, Rcyl, Nlcyl, Ncir, 1);
ver_cyl = ver_cyl(:,[3 1 2]);

Nvc = size(ver_cyl,1);


%%% indices of vertices at interface between cylinder and sphere
ind_int = (Nvc-Ncir+1):(Nvc);

ind_int = circshift(ind_int, +1*round(Nlcyl/2));

Iso = Nvc-Ncir;

%% create sphere section mesh

ver_sph = ver_cyl(ind_int, :);
fac_sph = zeros(0, 3);

tht_open = asin( Rcyl/Rsph );
tht_in   = pi-tht_open;
Larc = tht_in * Rsph;
Nlsph = ceil(Larc/delR_e);

delLs = Larc/Nlsph;
delTs = delLs/Rsph;

%%% latitude (polar) angles for sphere mesh
t1d = linspace(tht_open, pi, Nlsph+1);

Nc_cur = Ncir;
Nver_sph = 0;
p1d_cur = linspace(0, 2*pi, Ncir+1); p1d_cur = p1d_cur(1:end-1);

for itt = 1:(Nlsph)
    
    %%% check if need to halve number of facets
    Rnext = sin(t1d(itt+1))*Rsph;
    delC_next = 2*pi*Rnext/Nc_cur;
    
    if delC_next < RATIO_RC*delLs
        HalveFacets = 1;
    else
        HalveFacets = 0;
    end
    
    %%% create next ring vertices
    if HalveFacets==0
        Nc_next = Nc_cur;
        x1d = cos( t1d(itt+1) )*Rsph*ones(Nc_next,1);
        p1d_next = p1d_cur + (p1d_cur(2)-p1d_cur(1))/2;
        y1d = sin( t1d(itt+1) )*Rsph*cos( p1d_next );
        z1d = sin( t1d(itt+1) )*Rsph*sin( p1d_next );
    else
        Nc_next = Nc_cur/2;
        x1d = cos( t1d(itt+1) )*Rsph*ones(Nc_next,1);
        p1d_next = p1d_cur(1:2:end);
        y1d = sin( t1d(itt+1) )*Rsph*cos( p1d_next );
        z1d = sin( t1d(itt+1) )*Rsph*sin( p1d_next );
    end
    ver_sph = cat(1, ver_sph, [x1d(:) y1d(:) z1d(:)]);
    
    %%% create next ring faces
    if HalveFacets==0
        fac_tmp1 = Nver_sph +[ (1:Nc_cur)'                       , mod((1:Nc_cur)',Nc_cur)+1  , Nc_cur+(1:Nc_cur)'        ];
        fac_tmp2 = Nver_sph +[ mod((1:Nc_cur)',Nc_cur)+1+Nc_cur  , Nc_cur+(1:Nc_cur)'         , mod((1:Nc_cur)',Nc_cur)+1 ];
        fac_tmp = cat(1, fac_tmp1, fac_tmp2);
    else
        fac_tmp1 = Nver_sph + [ (1:2:Nc_cur)',    (2:2:Nc_cur)',                        Nc_cur+(1:Nc_next)'                ];
        fac_tmp2 = Nver_sph + [ (2:2:Nc_cur)',    1+mod((2:2:Nc_cur)',Nc_cur),          Nc_cur+1+mod((1:Nc_next)',Nc_next) ];
        fac_tmp3 = Nver_sph + [ (2:2:Nc_cur)',   Nc_cur+1+mod((1:Nc_next)',Nc_next),    Nc_cur+(1:Nc_next)'                ];
        fac_tmp = cat(1, fac_tmp1, fac_tmp2, fac_tmp3);
    end
    
    fac_sph = cat(1, fac_sph, fac_tmp);
    
    Nver_sph = Nver_sph + Nc_cur;
    
    Nc_cur = Nc_next;
    p1d_cur = p1d_next;
end
% figure(11); patchsexy('Faces',fac_sph,'Vertices',ver_sph);

TOL = 1e-6;
[ ver_fin, fac_fin ] = mergemeshes( ver_cyl, fac_cyl, ver_sph, fac_sph, TOL );
figure(12); patchsexy('Faces',fac_fin,'Vertices',ver_fin);

%% import mesh from stl file (generated in inventor)

[ fac_imp, ver_imp ] = stlread( './mesh_import/Helmet_v6p1_RFsurfonly_200119_v02_meter_bin.stl' );


figure(13); patchsexy('Faces',fac_imp,'Vertices',ver_imp);


%% interpolation onto shape of inventor mesh.


x_inv = ver_imp(:,1);
y_inv = ver_imp(:,2);
z_inv = ver_imp(:,3);

m_inv_cyl = ver_imp(:,1) >= -0.2;
m_inv_sph = ver_imp(:,1) <  0.2  & (x_inv<0.09346 | z_inv>-0.03914);

%%% cylinder points
x_inv_cyl = x_inv(m_inv_cyl);
y_inv_cyl = y_inv(m_inv_cyl);
z_inv_cyl = z_inv(m_inv_cyl);

ryz_inv_cyl = sqrt( y_inv_cyl.^2 + z_inv_cyl.^2 );
phi_inv_cyl = atan2( y_inv_cyl, z_inv_cyl );


%%% sphere points
x_inv_sph = x_inv(m_inv_sph);
y_inv_sph = y_inv(m_inv_sph);
z_inv_sph = z_inv(m_inv_sph);

ryz_inv_sph = sqrt( y_inv_sph.^2 + z_inv_sph.^2 );
r_inv_sph = sqrt( x_inv_sph.^2 + y_inv_sph.^2 + z_inv_sph.^2 );
phi_inv_sph = atan2( y_inv_sph, z_inv_sph );
tht_inv_sph = atan2( ryz_inv_sph, -x_inv_sph );

%%% make interp objects
SI_cyl = scatteredInterpolant( x_inv_cyl, y_inv_cyl, z_inv_cyl, ryz_inv_cyl,'linear','linear' );
% SI_cyl = scatteredInterpolant( x_inv_cyl, phi_inv_cyl, ryz_inv_cyl,'linear','none' );

% SI_sph = scatteredInterpolant( tht_inv_sph, y_inv_sph, z_inv_sph, r_inv_sph );
% SI_sph = scatteredInterpolant( tht_inv_sph, phi_inv_sph, r_inv_sph, 'linear','none');
SI_sph = scatteredInterpolant( tht_inv_sph.*cos(phi_inv_sph), tht_inv_sph.*sin(phi_inv_sph), r_inv_sph, 'linear','none');

%%%%%%%%%%%%%%%%%%%%%%%%% ab ovo mesh points, cyl and sph

x_fin = ver_fin(:,1);
y_fin = ver_fin(:,2);
z_fin = ver_fin(:,3);


m_fin_cyl = ver_fin(:,1) >= 0.2;
m_fin_sph = ver_fin(:,1) <  0.2 & (x_fin<0.09346 | z_fin>-0.03914);

%%% cylinder points
x_fin_cyl = x_fin(m_fin_cyl);
y_fin_cyl = y_fin(m_fin_cyl);
z_fin_cyl = z_fin(m_fin_cyl);



ryz_fin_cyl = sqrt( y_fin_cyl.^2 + z_fin_cyl.^2 );
phi_fin_cyl = atan2( y_fin_cyl, z_fin_cyl );

%%% sphere points
x_fin_sph = x_fin(m_fin_sph);
y_fin_sph = y_fin(m_fin_sph);
z_fin_sph = z_fin(m_fin_sph);

ryz_fin_sph = sqrt( y_fin_sph.^2 + z_fin_sph.^2 );
r_fin_sph = sqrt( x_fin_sph.^2 + y_fin_sph.^2 + z_fin_sph.^2 );
phi_fin_sph = atan2( y_fin_sph, z_fin_sph );
tht_fin_sph = atan2( ryz_fin_sph, -x_fin_sph );


r_fin_int_cyl = SI_cyl( x_fin_cyl, y_fin_cyl, z_fin_cyl );
% r_fin_int_cyl = SI_cyl( x_fin_cyl, phi_fin_cyl );

% r_fin_int_sph = SI_sph( tht_fin_sph, y_fin_sph, z_fin_sph );
% r_fin_int_sph = SI_sph( tht_fin_sph, phi_fin_sph );
r_fin_int_sph = SI_sph( tht_fin_sph.*cos(phi_fin_sph), tht_fin_sph.*sin(phi_fin_sph) );

%%%%%%%%%%%%%%%%%%%%%%%%% generate xyz coordinates for mesh

x_fin_int_cyl = x_fin_cyl;
y_fin_int_cyl = r_fin_int_cyl.* sin(phi_fin_cyl);
z_fin_int_cyl = r_fin_int_cyl.* cos(phi_fin_cyl);

x_fin_int_sph = -r_fin_int_sph .* cos(tht_fin_sph);
y_fin_int_sph =  r_fin_int_sph .* sin(tht_fin_sph) .* sin(phi_fin_sph);
z_fin_int_sph =  r_fin_int_sph .* sin(tht_fin_sph) .* cos(phi_fin_sph);

x_fin_int = nan(size(x_fin));
y_fin_int = nan(size(x_fin));
z_fin_int = nan(size(x_fin));

x_fin_int(m_fin_cyl) = x_fin_int_cyl; 
x_fin_int(m_fin_sph) = x_fin_int_sph;
y_fin_int(m_fin_cyl) = y_fin_int_cyl; 
y_fin_int(m_fin_sph) = y_fin_int_sph;
z_fin_int(m_fin_cyl) = z_fin_int_cyl; 
z_fin_int(m_fin_sph) = z_fin_int_sph;

ver_fin_int = zeros(size(ver_fin));
ver_fin_int(:,1) = x_fin_int;
ver_fin_int(:,2) = y_fin_int;
ver_fin_int(:,3) = z_fin_int;


figure(22); patchsexy('Faces',fac_fin,'Vertices',ver_fin_int);

% objmeshwrite(ver_fin_int, fac_fin, 'mesh/mesh_interp_200119_v01.obj');

[ ver_fin_nonan, fac_fin_nonan ] = mesh_removenanver( ver_fin_int, fac_fin );
figure(23); patchsexy('Faces',fac_fin_nonan,'Vertices',ver_fin_nonan);
objmeshwrite(ver_fin_nonan, fac_fin_nonan, './mesh/mesh_interp_200119_v01.obj');
