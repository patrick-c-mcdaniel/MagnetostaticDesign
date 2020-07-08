function [ Bx By Bz ] = SHcalc_L5( xx, yy, zz, m);

%     clear; close all;
%     load('precomp/MagROI_190220_v01.mat');
    
    i_m = m+6;


    xROI = xx(:);
    yROI = yy(:);
    zROI = zz(:);

    % N1d = 201;
    % 
    % x1d = linspace(-0.1,0.1,N1d);
    % z1d = linspace(-0.1,0.1,N1d);
    % 
    % [x2d z2d] = ndgrid(x1d, z1d);
    % 
    % xROI = x2d(:);
    % zROI = z2d(:);
    % yROI = zeros(size(xROI));
    % 
    xx_a = [ 0 1 ]';
    yy_a = [ 0 1 ]';
    zz_a = [ 0 1 ]';



    Nmag = numel(xx_a);
    Ntar = numel(xROI);

    MUZ = 4*pi*1e-7;   % vacuum permeability

    %% Create forward calculation dipole matrix
    %
    % [xROI yROI zROI] and [xx_a' yy_a' zz_a'] contain ROI and magnet
    %    coordinates, respectively. 
    % Create 

    D5_4d = zeros(3, 11, Ntar, Nmag);
    D5_2d = zeros(3*Ntar, 11*Nmag);

    % compute "outer difference" between sets of points:

    x_mt = repmat( xROI, [1 Nmag] ) - repmat( xx_a',[Ntar 1] );
    y_mt = repmat( yROI, [1 Nmag] ) - repmat( yy_a',[Ntar 1] );
    z_mt = repmat( zROI, [1 Nmag] ) - repmat( zz_a',[Ntar 1] );

    r_mt = sqrt( x_mt.^2 + y_mt.^2 + z_mt.^2 );   % r - target point distances (r) in magnet center coordinates
    t_mt = acos( z_mt ./ r_mt );                  % t - target point polar angles (theta) in magnet center coordinates
    p_mt = atan2( y_mt, x_mt );                   % p - target point azimuth angles (phi) in magnet center coordinates

    r = r_mt;
    rm7 = 1./(r.^7);
    ct = cos(t_mt);
    ct0 = ct.^0;
    ct1 = ct;
    ct2 = ct.^2;
    ct3 = ct.^3;
    ct4 = ct.^4;
    ct5 = ct.^5;
    ct6 = ct.^6;

    st = sin(t_mt);
    st0 = st.^0;
    st1 = st;
    st2 = st.^2;
    st3 = st.^3;
    st4 = st.^4;
    st5 = st.^5;
    st6 = st.^6;

    sp  = sin(1*p_mt);
    s1p = sin(1*p_mt);
    s2p = sin(2*p_mt);
    s3p = sin(3*p_mt);
    s4p = sin(4*p_mt);
    s5p = sin(5*p_mt);

    cp  = cos(1*p_mt);
    c1p = cos(1*p_mt);
    c2p = cos(2*p_mt);
    c3p = cos(3*p_mt);
    c4p = cos(4*p_mt);
    c5p = cos(5*p_mt);


    %%% rtp2xyz - coordinate transform matrices for (r,t,p)-hat basis vector
    %%% into (x,y,z)-hat basis vector



    % fill 4D L=5 multipole matrix (columns from m=-5 (index 1) to m=+5 (index 11)
    NORM = zeros(1,11);

    NORM(1)  = MUZ/(4*pi) * sqrt(4*pi/(11)) * 3/16*sqrt(77   /2/pi);
    NORM(2)  = MUZ/(4*pi) * sqrt(4*pi/(11)) * 3/16*sqrt(385  /pi  );
    NORM(3)  = MUZ/(4*pi) * sqrt(4*pi/(11)) * 1/16*sqrt(385  /2/pi);
    NORM(4)  = MUZ/(4*pi) * sqrt(4*pi/(11)) * 1/8 *sqrt(1155 /pi  );
    NORM(5)  = MUZ/(4*pi) * sqrt(4*pi/(11)) * 1/16*sqrt(165  /pi  );
    NORM(6)  = MUZ/(4*pi) * sqrt(4*pi/(11)) * 1/16*sqrt(11   /pi  );
    NORM(7)  = MUZ/(4*pi) * sqrt(4*pi/(11)) * 1/16*sqrt(165  /pi  );
    NORM(8)  = MUZ/(4*pi) * sqrt(4*pi/(11)) * 1/8 *sqrt(1155 /pi  );
    NORM(9)  = MUZ/(4*pi) * sqrt(4*pi/(11)) * 1/16*sqrt(385  /2/pi);
    NORM(10) = MUZ/(4*pi) * sqrt(4*pi/(11)) * 3/16*sqrt(385  /pi  );
    NORM(11) = MUZ/(4*pi) * sqrt(4*pi/(11)) * 3/16*sqrt(77   /2/pi);

    %% compute all 5th-order field terms using explicit analytical expressions
    %      That's right, I calculated all the 5th-order spherical harmonic
    %      field patterns by fucking hand.

    D5_4d(1, 1, :, :) = NORM(1)  * reshape( rm7.*(6*st6.*s5p.*cp                    - 5*st4.*ct2.*s5p.*cp                                                   + 5*st4.*c5p.*sp            ), [1 1 Ntar Nmag] );
    D5_4d(1, 2, :, :) = NORM(2)  * reshape( rm7.*(7*st5.*ct1.*s4p.*cp               - 4*st3.*ct3.*s4p.*cp                                                   + 4*st3.*ct1.*c4p.*sp       ), [1 1 Ntar Nmag] );
    D5_4d(1, 3, :, :) = NORM(3)  * reshape( rm7.*(6*st4.*(9*ct2-1).*s3p.*cp         - 3*st2.*ct2.*(9*ct2-1).*s3p.*cp    + 18*st4.*ct2.*s3p.*cp              + 3*st2.*(9*ct2-1).*c3p.*sp ), [1 1 Ntar Nmag] );
    D5_4d(1, 4, :, :) = NORM(4)  * reshape( rm7.*(6*st3.*(3*ct3-ct1).*s2p.*cp       - 2*st.*ct2.*(3*ct3-ct).*s2p.*cp    + st3.*ct.*(9*ct2-1).*s2p.*cp       + 2*st.*(3*ct3-ct).*c2p.*sp ), [1 1 Ntar Nmag] );
    D5_4d(1, 5, :, :) = NORM(5)  * reshape( rm7.*(6*st2.*(21*ct4-14*ct2+1).*sp.*cp  - ct2.*(21*ct4-14*ct2+1).*sp.*cp    + st2.*ct2.*(84*ct2-28).*sp.*cp     + (21*ct4-14*ct2+1).*cp.*sp ), [1 1 Ntar Nmag] );
    D5_4d(1, 6, :, :) = NORM(6)  * reshape( rm7.*(6*st.*(63*ct5-70*ct3+15*ct).*cp   + st.*ct.*(315*ct4-210*ct2+15).*cp                                                                  ), [1 1 Ntar Nmag] );
    D5_4d(1, 7, :, :) = NORM(7)  * reshape( rm7.*(6*st2.*(21*ct4-14*ct2+1).*cp.*cp  - ct2.*(21*ct4-14*ct2+1).*cp.*cp    + st2.*ct2.*(84*ct2-28).*cp.*cp     - (21*ct4-14*ct2+1).*sp.*sp ), [1 1 Ntar Nmag] );
    D5_4d(1, 8, :, :) = NORM(8)  * reshape( rm7.*(6*st3.*ct.*(3*ct2-1).*c2p.*cp     - 2*st.*ct3.*(3*ct2-1).*c2p.*cp     + st3.*ct.*(9*ct2-1).*c2p.*cp       - 2*st.*(3*ct3-ct).*s2p.*sp ), [1 1 Ntar Nmag] );
    D5_4d(1, 9, :, :) = NORM(9)  * reshape( rm7.*(6*st4.*(9*ct2-1).*c3p.*cp         - 3*st2.*ct2.*(9*ct2-1).*c3p.*cp    + 18*st4.*ct2.*c3p.*cp              - 3*st2.*(9*ct2-1).*s3p.*sp ), [1 1 Ntar Nmag] );
    D5_4d(1,10, :, :) = NORM(10) * reshape( rm7.*(7*st5.*ct.*c4p.*cp                - 4*st3.*ct3.*c4p.*cp                                                   - 4*st3.*ct.*s4p.*sp        ), [1 1 Ntar Nmag] );
    D5_4d(1,11, :, :) = NORM(11) * reshape( rm7.*(6*st6.*c5p.*cp                    - 5*st4.*ct2.*c5p.*cp                                                   - 5*st4.*s5p.*sp            ), [1 1 Ntar Nmag] );

    D5_4d(2, 1, :, :) = NORM(1)  * reshape( rm7.*(6*st6.*s5p.*sp                    - 5*st4.*ct2.*s5p.*sp                                                   - 5*st4.*c5p.*cp            ), [1 1 Ntar Nmag] );
    D5_4d(2, 2, :, :) = NORM(2)  * reshape( rm7.*(7*st5.*ct1.*s4p.*sp               - 4*st3.*ct3.*s4p.*sp                                                   - 4*st3.*ct1.*c4p.*cp       ), [1 1 Ntar Nmag] );
    D5_4d(2, 3, :, :) = NORM(3)  * reshape( rm7.*(6*st4.*(9*ct2-1).*s3p.*sp         - 3*st2.*ct2.*(9*ct2-1).*s3p.*sp    + 18*st4.*ct2.*s3p.*sp              - 3*st2.*(9*ct2-1).*c3p.*cp ), [1 1 Ntar Nmag] );
    D5_4d(2, 4, :, :) = NORM(4)  * reshape( rm7.*(6*st3.*(3*ct3-ct1).*s2p.*sp       - 2*st.*ct2.*(3*ct3-ct).*s2p.*sp    + st3.*ct.*(9*ct2-1).*s2p.*sp       - 2*st.*(3*ct3-ct).*c2p.*cp ), [1 1 Ntar Nmag] );
    D5_4d(2, 5, :, :) = NORM(5)  * reshape( rm7.*(6*st2.*(21*ct4-14*ct2+1).*sp.*sp  - ct2.*(21*ct4-14*ct2+1).*sp.*sp    + st2.*ct2.*(84*ct2-28).*sp.*sp     - (21*ct4-14*ct2+1).*cp.*cp ), [1 1 Ntar Nmag] );
    D5_4d(2, 6, :, :) = NORM(6)  * reshape( rm7.*(6*st.*(63*ct5-70*ct3+15*ct).*sp   + st.*ct.*(315*ct4-210*ct2+15).*sp                                                                  ), [1 1 Ntar Nmag] );
    D5_4d(2, 7, :, :) = NORM(7)  * reshape( rm7.*(6*st2.*(21*ct4-14*ct2+1).*cp.*sp  - ct2.*(21*ct4-14*ct2+1).*cp.*sp    + st2.*ct2.*(84*ct2-28).*cp.*sp     + (21*ct4-14*ct2+1).*sp.*cp ), [1 1 Ntar Nmag] );
    D5_4d(2, 8, :, :) = NORM(8)  * reshape( rm7.*(6*st3.*ct.*(3*ct2-1).*c2p.*sp     - 2*st.*ct3.*(3*ct2-1).*c2p.*sp     + st3.*ct.*(9*ct2-1).*c2p.*sp       + 2*st.*(3*ct3-ct).*s2p.*cp ), [1 1 Ntar Nmag] );
    D5_4d(2, 9, :, :) = NORM(9)  * reshape( rm7.*(6*st4.*(9*ct2-1).*c3p.*sp         - 3*st2.*ct2.*(9*ct2-1).*c3p.*sp    + 18*st4.*ct2.*c3p.*sp              + 3*st2.*(9*ct2-1).*s3p.*cp ), [1 1 Ntar Nmag] );
    D5_4d(2,10, :, :) = NORM(10) * reshape( rm7.*(7*st5.*ct.*c4p.*sp                - 4*st3.*ct3.*c4p.*sp                                                   + 4*st3.*ct.*s4p.*cp        ), [1 1 Ntar Nmag] );
    D5_4d(2,11, :, :) = NORM(11) * reshape( rm7.*(6*st6.*c5p.*sp                    - 5*st4.*ct2.*c5p.*sp                                                   + 5*st4.*s5p.*cp            ), [1 1 Ntar Nmag] );

    D5_4d(3, 1, :, :) = NORM(1)  * reshape( rm7.*(6*st5.*ct.*s5p                    + 5*st5.*ct.*s5p                                                        ), [1 1 Ntar Nmag] );
    D5_4d(3, 2, :, :) = NORM(2)  * reshape( rm7.*(6*st4.*ct2.*s4p                   - st6.*s4p                          + 4*st4.*ct2.*s4p                   ), [1 1 Ntar Nmag] );
    D5_4d(3, 3, :, :) = NORM(3)  * reshape( rm7.*(6*st3.*ct.*(9*ct2-1).*s3p         + 3*st3.*ct.*(9*ct2-1).*s3p         - 18*st5.*ct.*s3p                   ), [1 1 Ntar Nmag] );
    D5_4d(3, 4, :, :) = NORM(4)  * reshape( rm7.*(6*st2.*ct.*(3*ct3-ct).*s2p        + 2*st2.*ct.*(3*ct3-ct).*s2p        - st4.*(9*ct2-1).*s2p               ), [1 1 Ntar Nmag] );
    D5_4d(3, 5, :, :) = NORM(5)  * reshape( rm7.*(6*st.*ct.*(21*ct4-14*ct2+1).*sp   + st.*ct.*(21*ct4-14*ct2+1).*sp     - st3.*ct.*(84*ct2-28).*sp          ), [1 1 Ntar Nmag] );
    D5_4d(3, 6, :, :) = NORM(6)  * reshape( rm7.*(6*ct2.*(63*ct4-70*ct2+15)         - st2.*(315*ct4-210*ct2+15)                                             ), [1 1 Ntar Nmag] );
    D5_4d(3, 7, :, :) = NORM(7)  * reshape( rm7.*(6*st.*ct.*(21*ct4-14*ct2+1).*cp   + st.*ct.*(21*ct4-14*ct2+1).*cp     - st3.*ct.*(84*ct2-28).*cp          ), [1 1 Ntar Nmag] );
    D5_4d(3, 8, :, :) = NORM(8)  * reshape( rm7.*(6*st2.*ct.*(3*ct3-ct).*c2p        + 2*st2.*ct.*(3*ct3-ct).*c2p        - st4.*(9*ct2-1).*c2p               ), [1 1 Ntar Nmag] );
    D5_4d(3, 9, :, :) = NORM(9)  * reshape( rm7.*(6*st3.*ct.*(9*ct2-1).*c3p         + 3*st3.*ct.*(9*ct2-1).*c3p         - 18*st5.*ct.*c3p                   ), [1 1 Ntar Nmag] );
    D5_4d(3,10, :, :) = NORM(10) * reshape( rm7.*(6*st4.*ct2.*c4p                   - st6.*c4p                          + 4*st4.*ct2.*c4p                   ), [1 1 Ntar Nmag] );
    D5_4d(3,11, :, :) = NORM(11) * reshape( rm7.*(6*st5.*ct.*c5p                    + 5*st5.*ct.*c5p                                                        ), [1 1 Ntar Nmag] );

    % turn into 2D Dipole matrix (see slides to understand matrix structure)

    D5_2d(1:3:end, 1:11:end) = D5_4d(1, 1, :, :);
    D5_2d(1:3:end, 2:11:end) = D5_4d(1, 2, :, :);
    D5_2d(1:3:end, 3:11:end) = D5_4d(1, 3, :, :); 
    D5_2d(1:3:end, 4:11:end) = D5_4d(1, 4, :, :); 
    D5_2d(1:3:end, 5:11:end) = D5_4d(1, 5, :, :); 
    D5_2d(1:3:end, 6:11:end) = D5_4d(1, 6, :, :); 
    D5_2d(1:3:end, 7:11:end) = D5_4d(1, 7, :, :); 
    D5_2d(1:3:end, 8:11:end) = D5_4d(1, 8, :, :); 
    D5_2d(1:3:end, 9:11:end) = D5_4d(1, 9, :, :); 
    D5_2d(1:3:end,10:11:end) = D5_4d(1,10, :, :); 
    D5_2d(1:3:end,11:11:end) = D5_4d(1,11, :, :); 

    D5_2d(2:3:end, 1:11:end) = D5_4d(2, 1, :, :);
    D5_2d(2:3:end, 2:11:end) = D5_4d(2, 2, :, :);
    D5_2d(2:3:end, 3:11:end) = D5_4d(2, 3, :, :); 
    D5_2d(2:3:end, 4:11:end) = D5_4d(2, 4, :, :); 
    D5_2d(2:3:end, 5:11:end) = D5_4d(2, 5, :, :); 
    D5_2d(2:3:end, 6:11:end) = D5_4d(2, 6, :, :); 
    D5_2d(2:3:end, 7:11:end) = D5_4d(2, 7, :, :); 
    D5_2d(2:3:end, 8:11:end) = D5_4d(2, 8, :, :); 
    D5_2d(2:3:end, 9:11:end) = D5_4d(2, 9, :, :); 
    D5_2d(2:3:end,10:11:end) = D5_4d(2,10, :, :); 
    D5_2d(2:3:end,11:11:end) = D5_4d(2,11, :, :); 

    D5_2d(3:3:end, 1:11:end) = D5_4d(3, 1, :, :);
    D5_2d(3:3:end, 2:11:end) = D5_4d(3, 2, :, :);
    D5_2d(3:3:end, 3:11:end) = D5_4d(3, 3, :, :); 
    D5_2d(3:3:end, 4:11:end) = D5_4d(3, 4, :, :); 
    D5_2d(3:3:end, 5:11:end) = D5_4d(3, 5, :, :); 
    D5_2d(3:3:end, 6:11:end) = D5_4d(3, 6, :, :); 
    D5_2d(3:3:end, 7:11:end) = D5_4d(3, 7, :, :); 
    D5_2d(3:3:end, 8:11:end) = D5_4d(3, 8, :, :); 
    D5_2d(3:3:end, 9:11:end) = D5_4d(3, 9, :, :); 
    D5_2d(3:3:end,10:11:end) = D5_4d(3,10, :, :); 
    D5_2d(3:3:end,11:11:end) = D5_4d(3,11, :, :); 
    % 
    % save('D5matrix_190214_test.mat','D5_2d','D5_4d','Nmag','Ntar');

%     save('precomp/D5_190221_v03.mat','D5_2d','D5_4d','Nmag','Ntar');
    %% test...
    Mtest = zeros(22, 1);
    Mtest( i_m ) = 1;

    Bxyz = D5_2d * Mtest;
                
    Bx = Bxyz(1:3:end);
    By = Bxyz(2:3:end);
    Bz = Bxyz(3:3:end);
    
%     Bx_rs = reshape(Bx,[N1d, N1d]);
%     By_rs = reshape(By,[N1d, N1d]);
%     Bz_rs = reshape(Bz,[N1d, N1d]);
%     
%     figure(401); imagesc(Bx_rs); colormap jet; axis equal; caxis([-1e-4 1e-4]);
%     figure(402); imagesc(By_rs); colormap jet; axis equal; caxis([-1e-4 1e-4]);
%     figure(403); imagesc(Bz_rs); colormap jet; axis equal; caxis([-1e-4 1e-4]);
    %             
end