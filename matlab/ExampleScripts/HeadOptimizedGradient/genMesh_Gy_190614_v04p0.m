
Dcyl = 2*0.180;
Dsph = 2*0.195;
Dfpl = Dcyl - 2*0.03;

Rcyl = Dcyl/2;
Rsph = Dsph/2;
Rfpl = Dfpl/2;

Ncir = 2^6;

Ccyl = pi*Dcyl;
Csph = pi*Dsph;

delR_e = Ccyl/Ncir/2;

Xmax = 0.177;
Xbnd = ( Rsph^2 - Rcyl^2 )^(1/2);

Lcyl = Xmax-Xbnd;

Nlcyl = ceil(Lcyl/delR_e);

Nfpl = 8;

delL = Lcyl/Nlcyl;

RATIO_RC = 0.7;

%% create "front plate" mesh
% 
% [ ver_fpl, fac_fpl ] = genannumeshz( Xmax, Rfpl, Rcyl, Nfpl, Ncir, 1);
% ver_fpl = ver_fpl(:,[3 1 2]);

%% create cylinder mesh

[ ver_cyl, fac_cyl ] = gencylmeshz( Xmax, Xmax-Lcyl, Rcyl, Nlcyl, Ncir, 1);
ver_cyl = ver_cyl(:,[3 1 2]);

Nvc = size(ver_cyl,1);


%%% indices of vertices at interface between cylinder and sphere
ind_int = (Nvc-Ncir+1):(Nvc);
ind_int = circshift(ind_int, [0, round(Nlcyl/2)]);

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
% figure(10); patchsexy('Faces',fac_fpl,'Vertices',ver_fpl);

figure(11); patchsexy('Faces',fac_sph,'Vertices',ver_sph);

TOL = 1e-6;
[ ver_fin, fac_fin ] = mergemeshes( ver_cyl, fac_cyl, ver_sph, fac_sph, TOL );

% [ ver_fcl, fac_fcl ] = mergemeshes( ver_fin, fac_fin, ver_fpl, fac_fpl, TOL );

figure(12); patchsexy('Faces',fac_fin,'Vertices',ver_fin);
% figure(13); patchsexy('Faces',fac_fcl,'Vertices',ver_fcl);

objmeshwrite(ver_fin, fac_fin, './mesh/meshGy_Example_200626.obj');
