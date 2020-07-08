function [ output ] = proc_ROImask( fn_Mask, fn_optgeomsave )
    %% Notes: 190321
    %    - changed so that we only use the "shell" of the brain ROI for target
    %        points
    %% Notes: 190521
    %    - increase radii of magnet cylinder/sphere sections by 1cm
    %    - increase ROI extent towards foot by 1.4cm
    %% Notes: 190523
    %    - surface ROI mask only includes half of points (exploit symmetry for
    %        smaller computation
    %% Notes: 190524
    %   - added 30x head opening magnet locations (hopefully this will make it
    %       easier to get a homogeneous field down to the neck)
    %   - Turned script into function

    % mask_nii = load_nifti('/home3/patmcd/self_data/2017_10_19_Bay4/008/headmask_180808.nii.gz');
    % mask_vol = mask_nii.vol;
    % 
    % MSZ = size(mask_vol);
    % PSZ = [ 70 70 70 ];
    % 
    % [ x3 y3 z3 ] = ndgrid( 1:(2*PSZ(1)+MSZ(1)), 1:(2*PSZ(2)+MSZ(2)), 1:(2*PSZ(3)+MSZ(3) ));
    % 
    % mask_pad = zeros(size(mask_vol) + 2*PSZ);
    % mask_pad( (PSZ(1)+1):(MSZ(1)+PSZ(1)),  (PSZ(2)+1):(MSZ(2)+PSZ(2)),  (PSZ(3)+1):(MSZ(3)+PSZ(3)) ) = mask_vol;
    % 
    % MPSZ = size(mask_pad);
    % 
    % figure(3); imagesc( squeeze(mask_pad(floor(end/2),:,:)) ); colormap gray; axis image;
    % 
    % for islc = 1:MPSZ(1)
    %    
    %     mask_pad(islc,:,:) = reshape( imfill( squeeze( mask_pad(islc,:,:) ), 'holes' ), [1 MPSZ(2) MPSZ(3)] );
    %     
    %     
    % end
    % 
    % figure(13); imagesc( squeeze(mask_pad(floor(end/2),:,:)) ); colormap gray; axis image;
    % 
    % mask_dil = mask_pad;
    % 
    % Ndil = 35;
    % 
    % for idil = 1:1
    %     
    %     seball = strel('sphere',ceil(Ndil));
    %     mask_dil = imdilate( mask_dil, seball );
    %     
    % end
    % 
    % figure(21); imagesc( squeeze(mask_dil(:,:,floor(end/2))) ); colormap gray; axis image;
    % figure(22); imagesc( squeeze(mask_dil(:,floor(end/2),:)) ); colormap gray; axis image;
    % figure(23); imagesc( squeeze(mask_dil(floor(end/2),:,:)) ); colormap gray; axis image;

    % save('ROImask_180809.mat','mask_dil','mask_pad','MSZ','PSZ');

    %%% UNCOMMENT ABOVE TO LOAD/DILATE MASK

    load(fn_Mask);
    % load('../ROImask_180809.mat');
    [ x3 y3 z3 ] = ndgrid( 1:(2*PSZ(1)+MSZ(1)), 1:(2*PSZ(2)+MSZ(2)), 1:(2*PSZ(3)+MSZ(3) ));

    mask_dil2 = imdilate( mask_dil, ones([3 3 3]));
    mask_mag = mask_dil2 & (~mask_dil);
    mask_hd  = y3 <= MSZ(2)+PSZ(2);
    mask_mag_y = mask_mag & mask_hd;
    mask_mh  = (mask_mag_y | mask_pad);

    figure(33); imagesc( squeeze(mask_mh(floor(end/2),:,:)) ); colormap gray; axis image;

    %% mirror ROI, merge w/original, compute boundary of joint

    mask_dil_mir = flip( mask_dil, 3);
    mask_dil_or  = mask_dil | mask_dil_mir;
    mask_dil_or2 = imdilate( mask_dil_or, ones([3 3 3]));

    mask_mag_or  = mask_dil_or2 & (~mask_dil_or);

    mask_mag_or_y = mask_mag_or & mask_hd;
    mask_moy_h   = mask_mag_or_y | mask_pad;

    figure(43); imagesc( squeeze(mask_moy_h(floor(end/2),:,:)) ); colormap gray; axis image;


    % save('mask_mag_or_y.mat','mask_mag_or_y');
    close all
    %% look at magnet surface


    figure(53); imagesc( squeeze(mask_mag_or_y(floor(end/2),:,:)) ); colormap gray; axis image;

    indpos = find( mask_mag_or_y );
    [xp zp yp] = ind2sub( size(mask_mag_or_y), indpos );

    %%% define center of head/origin for coordinate conversion
    xyzor = [ 187 158.5 158.5 ];
    [ lp pp rp RMAX ] = xyz_to_lpr( xp, yp, zp, xyzor );

    figure(55); scatter3( xp, yp, zp, 1, lp ); axis equal; colormap jet; colorbar;
    figure(56); scatter3( xp, yp, zp, 1, pp); axis equal; colormap jet; colorbar;
    figure(57); scatter3( xp, yp, zp, 1, rp); axis equal; colormap jet; colorbar;
    close( [53 55 56 57] );

    % figure(56); plot3( rp.*cos(pp), rp.*sin(pp), lp,'.'); axis equal
    %% mirror head mask to get target field ROI

    mask_pad_mir = flip( mask_pad, 3 );
    mask_pad_or  = mask_pad | mask_pad_mir;

    mask_pad_or_sm = imerode( imdilate( mask_pad_or, strel('sphere',4) ), strel('sphere',4) );

    figure(61); imagesc( squeeze(mask_pad_or_sm(:,:,floor(end/2))) ); colormap gray; axis image;
    figure(62); imagesc( squeeze(mask_pad_or_sm(:,floor(end/2),:)) ); colormap gray; axis image;
    figure(63); imagesc( squeeze(mask_pad_or_sm(floor(end/2),:,:)) ); colormap gray; axis image;

    %% define magnet center points (in "physical" coordinates)
    %  
    %   notes on global coordinate system used:
    %       - all values in [m] (METERS)
    %       - origin is at center of hemisphere == center of "head end" of
    %           cylinder
    %       - +z axis: transverse to cylinder, towards "top" of sphere, along
    %           +B0
    %       - +x axis: along cylinder axis, pointing towards patient
    %       - +y axis: defined by +x and +z axes; points transverse to cylinder
    %           perpendicular to B0
    %       - This is NOT the coordinate system used in the mask files, since
    %           that coordinate system was inherited from the base *.nii volumes,
    %           which were in turn produced on a typical MRI system
    %
    %   notes on cylinder coordinate system used:
    %       - for convenience, use a 

    %% Cylinder points + initial guess magnetization
    MUZ    = 4*pi*1e-7;                         % vacuum permittivity (SI units)
    Mcyl   = 1.45 * (0.0254^2)*0.0127 / MUZ;           %  Br * Vblk / MUZ = block magnetic dipole moment magnitude

    Nc_tht = 26;            % number of blocks in single ring

    % Rcyl   = 0.302/2;       % meters
    % Rcyl2  = Rcyl + 0.015;   % bigger bulb radius [m]

    %%% edit 19/05/21
    Rcyl_sm = 0.302/2;
    Rcyl   = 0.302/2 + 0.01; % meters
    Rcyl2  = Rcyl + 0.015;   % meters

    Lcyl   = 0.1675;    
    Nc_len = 7;             

    ll1d   = linspace(0, 1*Lcyl, Nc_len);                                      % set of l-coordinates of cylinder rings
    tt1d   = linspace(0, 2*pi, Nc_tht+1); tt1d = tt1d(1:(end-1));               % set of theta-coordinates of cylinder rungs

    [ ll2d tt2d ] = ndgrid(ll1d, tt1d);

    %%% cylinder radius varies as we get near the "bulb"

    ll_cut = sqrt(Rcyl2^2 - Rcyl^2);

    Rcyl_l = Rcyl * (ll2d>=ll_cut) + sqrt((Rcyl2^2 - ll2d.^2)) .* (ll2d<ll_cut);

    xx_c = ll2d;
    yy_c = Rcyl_l .* sin(tt2d);
    zz_c = Rcyl_l .* cos(tt2d);

    %%%% cylinder block orientations (only for initial guess)
    % 
    % for magnet orientation, use normal spherical coordinates:
    %       theta = spherical angle off of +z
    %       phi   = azimuthal angle off +x, increasing towards +y

    tt_bc = 2*tt2d;                         % Uniform halbach magnetization rotates twice as you move around ring
    pp_bc = pi/2 * ones(size(tt2d));        % no x-component of any magnetization
    mm_bc = Mcyl * ones(size(tt2d));

    %%%%%%% extra cylinder end points

    x1d_op = 0.1524 - 0.0254*[0 1 2];

    y1d_op = [  0.117,  0.132,  0.132,  0.132,  0.117 , ...
               -0.117, -0.132, -0.132, -0.132, -0.117 ];

    z1d_op = [  0.068,  0.037,  0.000, -0.037, -0.068 , ...
                0.068,  0.037,  0.000, -0.037, -0.068 ];

    [ xc_op, yc_op ] = ndgrid( x1d_op, y1d_op );
    [ xc_op, zc_op ] = ndgrid( x1d_op, z1d_op );

    tc_op = zeros(size(xc_op));
    pc_op = zeros(size(xc_op));
    mc_op = zeros(size(xc_op));

    %% sphere points + inital guess magnetization
    Msph = 3/4*Mcyl;                        % sphere has 4/3 field of cylinder, so scale down magnetization by 3/4 for uniform field

    Ns_tht = floor( (Nc_tht-1)/2 );         % number of "rungs" in spherical section (not counting theta=0 or theta=pi)
    DELtht = 2*pi/Nc_tht;
    DELblk = 0.032;                         % block center approximate spacing along arc (appx==DELphi*Rs1d(theta))

    ts1d   = DELtht:DELtht:(Ns_tht*DELtht);
    Rsph   = Rcyl;
    Rsph2  = Rcyl2;
    % Rs1d   = Rsph * sin(ts1d);
    % C2s1d  = pi * Rs1d;
    % Ns_phs = ceil( C2s1d / DELblk ) - 1;            % array with Ns_tht entries, each with the number of phi-values at that theta-index

    %%% edit 19/05/21
    Rs1d   = Rcyl_sm * sin(ts1d);
    C2s1d  = pi * Rs1d;
    Ns_phs = ceil( C2s1d / DELblk ) - 1;            % array with Ns_tht entries, each with the number of phi-values at that theta-index

    DELphi = pi ./ (Ns_phs+1);                      % array of phi step sizes at each theta-level

    tt_s = zeros(1,0);
    pp_s = zeros(1,0);

    for iss = 1:Ns_tht
        ps1d = DELphi(iss):DELphi(iss):(Ns_phs(iss)*DELphi(iss));

        for ipp = 1:Ns_phs(iss)
            tt_s = cat(2, tt_s, ts1d(iss));
            pp_s = cat(2, pp_s, ps1d(ipp));         % note: for pp_s, phi is defined as the angle from +y, initially increasing towards -x
                                                    % it's peculiar, but it is the most natural choice for this part of the magnet


        end
    end

    xx_s = Rsph2 * sin(tt_s) .* sin(-1*pp_s);
    yy_s = Rsph2 * sin(tt_s) .* cos(-1*pp_s);
    zz_s = Rsph2 * cos(tt_s);

    %%% Get all (theta, phi) coordinates for blocks (only for initial guess)

    tt_bs = 2*tt_s;                         % spherical halbach has same theta-dependence as cylindrical
    pp_bs = pi/2 + pp_s;                    % same as block position, but need to add pi/2 to convert to global spherical coordinates
    mm_bs = Msph * ones(size(tt_s));

    %% get all (x,y,z) coordinates, and inital block (M, theta, phi) coordinates together

    xx_a = cat(1, xx_c(:), xx_s(:), xc_op(:));
    yy_a = cat(1, yy_c(:), yy_s(:), yc_op(:));
    zz_a = cat(1, zz_c(:), zz_s(:), zc_op(:));

    tt_ba = cat(1, tt_bc(:), tt_bs(:), tc_op(:));
    pp_ba = cat(1, pp_bc(:), pp_bs(:), pc_op(:));
    mm_ba = cat(1, mm_bc(:), mm_bs(:), mc_op(:));

    mx_ba = mm_ba .* sin(tt_ba) .* cos(pp_ba);
    my_ba = mm_ba .* sin(tt_ba) .* sin(pp_ba);
    mz_ba = mm_ba .* cos(tt_ba);

    figure(101); hold off; quiver3( xx_a, yy_a, zz_a, mx_ba, my_ba, mz_ba ); axis equal

    %% figure out target point coordinates
    %
    %   These are points inside of the head where we want to optimize the
    %   magnetic field. These coordinates will include the brain (where we
    %   ultimately wish to image), skull, and scalp, and go down to the
    %   axial->cor plane that goes through the medulla/spine interface and is
    %   tangent to the bottom of the eyeballs (arbitrary definitions)

    DELtar = 7e-3;         % target point spacing (isotropic; units [m])

    xy1 = [ 118 200 0 ];                                    % points defining limit plane for target ROI
    xy2 = [ 212 228 0 ];                                    % correspont to lower eye + lower medulla
    v12_n = [ 0  1  0; -1  0  0; 0  0  1 ] * (xy1'-xy2');   % vector normal to limit plane
    % v12_n = [-28 94 0];

    x1tar = 0:DELtar:0.15; x1tar = cat(2, -1*flip(x1tar(2:end)), x1tar);
    y1tar = 0:DELtar:0.15; y1tar = cat(2, -1*flip(y1tar(2:end)), y1tar);
    z1tar = 0:DELtar:0.15; z1tar = cat(2, -1*flip(z1tar(2:end)), z1tar);

    [x3tar y3tar z3tar] = ndgrid(x1tar, y1tar, z1tar);

    %%% convert to head nii units
    % xyzor = [ 187 158.5 132 ]
    %
    x3nii = z3tar*1000 + xyzor(1);
    y3nii = x3tar*1000 + xyzor(2);
    z3nii = y3tar*1000 + xyzor(3);

    ROImsk = (x3nii - xy1(1))*v12_n(1) + (y3nii - xy1(2))*v12_n(2) <= 0;

    x3ind = max(min( size(mask_pad_or_sm,1), round(x3nii)), 1);
    y3ind = max(min( size(mask_pad_or_sm,1), round(y3nii)), 1);
    z3ind = max(min( size(mask_pad_or_sm,1), round(z3nii)), 1);

    niimsk = reshape( mask_pad_or_sm( sub2ind( size(mask_pad_or_sm), x3ind(:), y3ind(:), z3ind(:) ) )==0, size(x3ind));

    ROImsk( niimsk ) = 0;

    %%% "slide down" to give extra ~1.5cm to ~2cm of ROI; in case subject's
    %%% neck is short.
    %%%   Note: hard-coded 2x mask shifts gives 1.4cm extra; different ROI
    %%%   resolution will require different number of shifts

    ROImsk_s1 = circshift(ROImsk,[1 0 0]);
    ROImsk_s2 = circshift(ROImsk,[2 0 0]);

    ROImsk = ROImsk | ROImsk_s1 | ROImsk_s2;

    %%% ROImsk_or defines actual target ROI
    %     (x3tar, y3tar, z3tar) are coorindates of potential target points
    %     [x3tar, y3tar, z3tar](ROImsk_or==1) == actual target point
    %           coordinates
    %
    ROImsk_or = ROImsk | flip(ROImsk,2);
    %% generate "shell" mask from 3D volume mask
    ROImsk_or_vol = ROImsk_or;
    ROImsk_or_ero = imerode(ROImsk_or_vol,ones([3 3 3]));
    ROImsk_or = ROImsk_or_vol & (~ROImsk_or_ero);

    %%% only look at half of ROImsk_or boundary:

    ROImsk_or  = ROImsk_or & y3tar>=0;

    figure(101); hold on;
    plot3(x3tar(ROImsk_or), y3tar(ROImsk_or), z3tar(ROImsk_or),'.')
    %%%

    %% block indices for midline, ypos, yneg regions
    inds_all = 1:numel(mx_ba);
    inds_mid = [ 1:7 184 196 219 260 283 295 92:98 ];

    %%% cylinder mirror indices
    inds_pos_cyl = [    8:14    , ...
                       15:21    , ...
                       22:28    , ...
                       29:35    , ...
                       36:42    , ...
                       43:49    , ...
                       50:56    , ...
                       57:63    , ...
                       64:70    , ...
                       71:77    , ...
                       78:84    , ...
                       85:91    , ...
                      297:299   , ...   %%% add new end ring blocks
                      300:302   , ...
                      303:305   , ...
                      306:308   , ...
                      309:311   ];

    inds_neg_cyl = [  176:182   , ...
                      169:175   , ...
                      162:168   , ...
                      155:161   , ...
                      148:154   , ...
                      141:147   , ...
                      134:140   , ...
                      127:133   , ...
                      120:126   , ...
                      113:119   , ...
                      106:112   , ...
                       99:105   , ...
                      312:314   , ...   %%% add new end ring blocks
                      315:317   , ...
                      318:320   , ...
                      321:323   , ...
                      324:326   ];

    %%% sphere mirror indices
    inds_pos_sph = [  183           , ...
                      188:-1:186    , ...
                      195:-1:192    , ...
                      206:-1:201    , ...
                      218:-1:213    , ...
                      232:-1:226    , ...
                      246:-1:240    , ...
                      259:-1:254    , ...
                      272:-1:267    , ...
                      282:-1:279    , ...
                      290:-1:288    , ...
                      294           ];
    inds_neg_sph = [  185           , ...
                      189:1:191    , ...
                      197:1:200    , ...
                      207:1:212    , ...
                      220:1:225    , ...
                      233:1:239    , ...
                      247:1:253    , ...
                      261:1:266    , ...
                      273:1:278    , ...
                      284:1:287    , ...
                      291:1:293    , ...
                      296           ];

    %%% unite cylinder+sphere indices
    inds_pos = [ inds_pos_cyl inds_pos_sph ];
    inds_neg = [ inds_neg_cyl inds_neg_sph ];

    % dispimgs(double(ROImsk_or),'gray','axial');
    save(fn_optgeomsave, ...                                                       % filename
         'xx_a','yy_a','zz_a','tt_ba','pp_ba','mm_ba','mx_ba','my_ba','mz_ba', ...      % parameters to define magnet positions/initial guess
         'ROImsk_or','ROImsk_or_vol','x3tar','y3tar','z3tar', ...                       % ROI coordinate matrices + ROI binary mask
         'inds_mid','inds_pos','inds_neg');                                          
     output = 1;
end