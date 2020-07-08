%% process 200119_v05 (y-pol, w/zero ROI, lambda=6.0e-6) into contours

% str_200119_v05 = 'RF_y_zero_r08p5_200119_v05';

fn_coil_red = './optresults/coil_red_RFLR_Example_200626.mat';
    
BSFN = [ './BSfiles/BS_RFLR_Example_200626.biot' ];
 
load(fn_coil_red);

%%%

displayStreamFunction(coil_red.listTriangle,coil_red.s,coil_red.listNode)

ver_mesh = coil_red.listNode;
fac_mesh = coil_red.listTriangle;
psi_mesh = coil_red.s;
xyz_ori = [ 0.1 0 0 ];
tz_bas = [ 0 0 -1 ; 0 -1 0 ; -1 0 0 ];
Niso = 16;


[ cont_xyz ] = genisocont_sphcord_v2( ver_mesh, fac_mesh, psi_mesh, xyz_ori, tz_bas, Niso, 'BSFN', BSFN, 'PLOT', 1, 'stdSM', 35 );

%% remove spurious contours

Ncont = numel(cont_xyz);

figure(22); hold off;

for icc = [ 1 8 9 11 12 13 15 22 ]
    cont_tmp = cont_xyz{icc};
    plot3( cont_tmp(1,:), cont_tmp(2,:), cont_tmp(3,:) ); grid on; grid minor; axis equal;
    hold on;
end

cont_xyz_rem = cell(Ncont - 8, 1);
cc_use = [2:7 10 14 16:21 ];

for icc = 1:numel(cont_xyz_rem)
   cont_xyz_rem{icc} = cont_xyz{cc_use(icc)};
    
end

writeBScont_xyz(cont_xyz_rem, BSFN);    

%% symmetrize contours - use yneg contours; mirror

figure(31); hold off;

for icc = [ 1:7 ]
    cont_tmp = cont_xyz_rem{icc};
    plot3( cont_tmp(1,:), cont_tmp(2,:), cont_tmp(3,:) ); grid on; grid minor; axis equal;
    hold on;
    cont_xyz_rem{numel(cont_xyz_rem)-icc+1} = cat(1, cont_tmp(1,:), -cont_tmp(2,:), cont_tmp(3,:));
end

writeBScont_xyz(cont_xyz_rem, BSFN);    
save('./optresults/cont_xyz_RFLR_Example_200626.mat','cont_xyz_rem');

%% load nose surface and use for generating z-offset near nose ridge

[ fac_nose, ver_nose ] = stlread( './STL_import/Helmet_v6p1_nose_v01_200120_meter_bin.stl' );

SI_nose = scatteredInterpolant( ver_nose(:, 1), ver_nose(:,2), ver_nose(:,3) );

figure(41); hold off;
       

cont_xyz_nose = cell(size(cont_xyz_rem));

for icc = 1:14
    cont_tmp = cont_xyz_rem{icc};
    x_tmp = cont_tmp(1,:);
    y_tmp = cont_tmp(2,:);
    z_tmp = cont_tmp(3,:);
    msk_nose = (z_tmp < 0) & ...
               ( ( ( x_tmp > 1.214e-3 ) & ( (x_tmp - 1.214e-3) / (50.744e-3) - (abs(y_tmp) / (30.6e-3) ) > 0 ) ) & (x_tmp <= 51.958e-3) | ...
                 ( ( x_tmp > 51.958e-3) & ( abs(y_tmp) < 30.6e-3 ) ) );
       
    z_nose = SI_nose(  x_tmp(msk_nose), y_tmp(msk_nose) );
    z_tmp(msk_nose) = z_nose + 0.7e-3; 
    cont_tmp_nose = cat(1, x_tmp, y_tmp, z_tmp );
    
    plot3( cont_tmp_nose(1,:), cont_tmp_nose(2,:), cont_tmp_nose(3,:) ); grid on; grid minor; axis equal;
    hold on;
    
    cont_xyz_nose{icc} = cont_tmp_nose;
    
end
save('./optresults/cont_xyz_nose_RFLR_Example_200626.mat','cont_xyz_nose');
%% combine RO/shim coil with magnet model in single BS-file
% 
% BSFN_mag = '/home/patmcd/Documents/HelmetMagnet/Magnet_v2/190327_proc/BSfiles/Magnet_190327_v01.biot';
% 
% BSFN_mag_ShimRO = 'BSfiles/Mag_ShimROy_190402_v04_test.biot';
% 
% BSfile_combine( BSFN_mag, BSFN, BSFN_mag_ShimRO );
