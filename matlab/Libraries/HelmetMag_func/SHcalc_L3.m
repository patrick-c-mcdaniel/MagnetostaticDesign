function [ Bx By Bz ] = SHcalc_L3( xx, yy, zz, m);

%     clear; close all;
%     load('precomp/MagROI_190220_v01.mat');
    
    i_m = m+4;


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

    D3_4d = zeros(3, 7, Ntar, Nmag);
    D3_2d = zeros(3*Ntar, 7*Nmag);

    % compute "outer difference" between sets of points:

    x_mt = repmat( xROI, [1 Nmag] ) - repmat( xx_a',[Ntar 1] );
    y_mt = repmat( yROI, [1 Nmag] ) - repmat( yy_a',[Ntar 1] );
    z_mt = repmat( zROI, [1 Nmag] ) - repmat( zz_a',[Ntar 1] );

    r_mt = sqrt( x_mt.^2 + y_mt.^2 + z_mt.^2 );   % r - target point distances (r) in magnet center coordinates
    t_mt = acos( z_mt ./ r_mt );                  % t - target point polar angles (theta) in magnet center coordinates
    p_mt = atan2( y_mt, x_mt );                   % p - target point azimuth angles (phi) in magnet center coordinates

    r = r_mt;
    rm5 = 1./(r.^5);
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
    NORM = zeros(1,7);

    % NORM(1)  = MUZ/(4*pi) * sqrt(4*pi/(7)) * 1/4 *sqrt(35   /2/pi);
    % NORM(2)  = MUZ/(4*pi) * sqrt(4*pi/(7)) * 1/4 *sqrt(105  /pi  );
    % NORM(3)  = MUZ/(4*pi) * sqrt(4*pi/(7)) * 1/4 *sqrt(21   /2/pi);
    % NORM(4)  = MUZ/(4*pi) * sqrt(4*pi/(7)) * 1/4 *sqrt(7    /pi  );
    % NORM(5)  = MUZ/(4*pi) * sqrt(4*pi/(7)) * 1/4 *sqrt(21   /2/pi);
    % NORM(6)  = MUZ/(4*pi) * sqrt(4*pi/(7)) * 1/4 *sqrt(105  /pi  );
    % NORM(7)  = MUZ/(4*pi) * sqrt(4*pi/(7)) * 1/4 *sqrt(35   /2/pi);

    NORM(1)  = MUZ/(4*pi) * sqrt(4*pi/(7)) * 1/4 *sqrt(35   /2/pi);
    NORM(2)  = MUZ/(4*pi) * sqrt(4*pi/(7)) * 1/4 *sqrt(105  /pi  );
    NORM(3)  = MUZ/(4*pi) * sqrt(4*pi/(7)) * 1/4 *sqrt(21   /2/pi);
    NORM(4)  = MUZ/(4*pi) * sqrt(4*pi/(7)) * 1/4 *sqrt(7    /pi  );
    NORM(5)  = MUZ/(4*pi) * sqrt(4*pi/(7)) * 1/4 *sqrt(21   /2/pi);
    NORM(6)  = MUZ/(4*pi) * sqrt(4*pi/(7)) * 1/4 *sqrt(105  /pi  );
    NORM(7)  = MUZ/(4*pi) * sqrt(4*pi/(7)) * 1/4 *sqrt(35   /2/pi);

    %% compute all 3rd-order field terms using explicit analytical expressions

    D3_4d(1, 1, :, :) = NORM(1)  * reshape( rm5.*(4*st4.*s3p.*cp                    - 3*st2.*ct2.*s3p.*cp                                                   + 3*st2.*c3p.*sp            ), [1 1 Ntar Nmag] );
    D3_4d(1, 2, :, :) = NORM(2)  * reshape( rm5.*(4*st3.*ct.*s2p.*cp                - 2*st.*ct3.*s2p.*cp                + st3.*ct.*s2p.*cp                  + 2*st.*ct.*c2p.*sp         ), [1 1 Ntar Nmag] );
    D3_4d(1, 3, :, :) = NORM(3)  * reshape( rm5.*((16*st2.*ct2-4*st4).*sp.*cp       - (4*ct4-11*st2.*ct2).*sp.*cp                                           + (4*ct2-st2).*cp.*sp       ), [1 1 Ntar Nmag] );
    D3_4d(1, 4, :, :) = NORM(4)  * reshape( rm5.*(st.*(20*ct3-12*ct).*cp            + ct.*(15*st.*ct2-3*st).*cp                                                                         ), [1 1 Ntar Nmag] );
    D3_4d(1, 5, :, :) = NORM(5)  * reshape( rm5.*((16*st2.*ct2-4*st4).*cp.*cp       - (4*ct4-11*st2.*ct2).*cp.*cp                                           - (4*ct2-st2).*sp.*sp       ), [1 1 Ntar Nmag] );
    D3_4d(1, 6, :, :) = NORM(6)  * reshape( rm5.*(4*st3.*ct.*c2p.*cp                - 2*st.*ct3.*c2p.*cp                + st3.*ct.*c2p.*cp                  - 2*st.*ct.*s2p.*sp         ), [1 1 Ntar Nmag] );
    D3_4d(1, 7, :, :) = NORM(7)  * reshape( rm5.*(4*st4.*c3p.*cp                    - 3*st2.*ct2.*c3p.*cp                                                   - 3*st2.*s3p.*sp            ), [1 1 Ntar Nmag] );

    D3_4d(2, 1, :, :) = NORM(1)  * reshape( rm5.*(4*st4.*s3p.*sp                    - 3*st2.*ct2.*s3p.*sp                                                   - 3*st2.*c3p.*cp            ), [1 1 Ntar Nmag] );
    D3_4d(2, 2, :, :) = NORM(2)  * reshape( rm5.*(4*st3.*ct.*s2p.*sp                - 2*st.*ct3.*s2p.*sp                + st3.*ct.*s2p.*sp                  - 2*st.*ct.*c2p.*cp         ), [1 1 Ntar Nmag] );
    D3_4d(2, 3, :, :) = NORM(3)  * reshape( rm5.*((16*st2.*ct2-4*st4).*sp.*sp       - (4*ct4-11*st2.*ct2).*sp.*sp                                           - (4*ct2-st2).*cp.*cp       ), [1 1 Ntar Nmag] );
    D3_4d(2, 4, :, :) = NORM(4)  * reshape( rm5.*(st.*(20*ct3-12*ct).*sp            + ct.*(15*st.*ct2-3*st).*sp                                                                         ), [1 1 Ntar Nmag] );
    D3_4d(2, 5, :, :) = NORM(5)  * reshape( rm5.*((16*st2.*ct2-4*st4).*cp.*sp       - (4*ct4-11*st2.*ct2).*cp.*sp                                           + (4*ct2-st2).*sp.*cp       ), [1 1 Ntar Nmag] );
    D3_4d(2, 6, :, :) = NORM(6)  * reshape( rm5.*(4*st3.*ct.*c2p.*sp                - 2*st.*ct3.*c2p.*sp                + st3.*ct.*c2p.*sp                  + 2*st.*ct.*s2p.*cp         ), [1 1 Ntar Nmag] );
    D3_4d(2, 7, :, :) = NORM(7)  * reshape( rm5.*(4*st4.*c3p.*sp                    - 3*st2.*ct2.*c3p.*sp                                                   + 3*st2.*s3p.*cp            ), [1 1 Ntar Nmag] );

    D3_4d(3, 1, :, :) = NORM(1)  * reshape( rm5.*(4*st3.*ct.*s3p                    + 3*st3.*ct.*s3p                                                        ), [1 1 Ntar Nmag] );
    D3_4d(3, 2, :, :) = NORM(2)  * reshape( rm5.*(4*st2.*ct2.*s2p                   - st4.*s2p                          + 2*st2.*ct2.*s2p                   ), [1 1 Ntar Nmag] );
    D3_4d(3, 3, :, :) = NORM(3)  * reshape( rm5.*((16*st.*ct3-4*st3.*ct).*sp        + (4*st.*ct3-11*st3.*ct).*sp                                            ), [1 1 Ntar Nmag] );
    D3_4d(3, 4, :, :) = NORM(4)  * reshape( rm5.*(ct.*(20*ct3-12*ct)                - st.*(15*st.*ct2-3*st)                                                 ), [1 1 Ntar Nmag] );
    D3_4d(3, 5, :, :) = NORM(5)  * reshape( rm5.*((16*st.*ct3-4*st3.*ct).*cp        + (4*st.*ct3-11*st3.*ct).*cp                                            ), [1 1 Ntar Nmag] );
    D3_4d(3, 6, :, :) = NORM(6)  * reshape( rm5.*(4*st2.*ct2.*c2p                   - st4.*c2p                          + 2*st2.*ct2.*c2p                   ), [1 1 Ntar Nmag] );
    D3_4d(3, 7, :, :) = NORM(7)  * reshape( rm5.*(4*st3.*ct.*c3p                    + 3*st3.*ct.*c3p                                                        ), [1 1 Ntar Nmag] );

    % turn into 2D Dipole matrix (see slides to understand matrix structure)

    D3_2d(1:3:end, 1:7:end) = D3_4d(1, 1, :, :);
    D3_2d(1:3:end, 2:7:end) = D3_4d(1, 2, :, :);
    D3_2d(1:3:end, 3:7:end) = D3_4d(1, 3, :, :); 
    D3_2d(1:3:end, 4:7:end) = D3_4d(1, 4, :, :); 
    D3_2d(1:3:end, 5:7:end) = D3_4d(1, 5, :, :); 
    D3_2d(1:3:end, 6:7:end) = D3_4d(1, 6, :, :); 
    D3_2d(1:3:end, 7:7:end) = D3_4d(1, 7, :, :); 

    D3_2d(2:3:end, 1:7:end) = D3_4d(2, 1, :, :);
    D3_2d(2:3:end, 2:7:end) = D3_4d(2, 2, :, :);
    D3_2d(2:3:end, 3:7:end) = D3_4d(2, 3, :, :); 
    D3_2d(2:3:end, 4:7:end) = D3_4d(2, 4, :, :); 
    D3_2d(2:3:end, 5:7:end) = D3_4d(2, 5, :, :); 
    D3_2d(2:3:end, 6:7:end) = D3_4d(2, 6, :, :); 
    D3_2d(2:3:end, 7:7:end) = D3_4d(2, 7, :, :); 

    D3_2d(3:3:end, 1:7:end) = D3_4d(3, 1, :, :);
    D3_2d(3:3:end, 2:7:end) = D3_4d(3, 2, :, :);
    D3_2d(3:3:end, 3:7:end) = D3_4d(3, 3, :, :); 
    D3_2d(3:3:end, 4:7:end) = D3_4d(3, 4, :, :); 
    D3_2d(3:3:end, 5:7:end) = D3_4d(3, 5, :, :); 
    D3_2d(3:3:end, 6:7:end) = D3_4d(3, 6, :, :); 
    D3_2d(3:3:end, 7:7:end) = D3_4d(3, 7, :, :); 

    %% test...
    Mtest = zeros(14, 1);
    Mtest( i_m ) = 1;

    Bxyz = D3_2d * Mtest;
                
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