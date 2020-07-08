
load('./optresults/Gx_coil_red_Example_200626.mat');

%% interp. vertices (at native mesh resolution)

[ ver_up, fac_up, psi_up ] = meshScalarUpsample( coil_red.listNode, coil_red.listTriangle, coil_red.s, 3 );

%% generate isocontours


displayStreamFunction(fac_up,psi_up,ver_up);

ver_mesh = ver_up;
fac_mesh = fac_up;
psi_mesh = psi_up;
xyz_ori = [ 0.15 0 0 ];
tz_bas = [ 0 0 -1 ; 0 -1 0 ; -1 0 0 ];
Niso = 30;

BSFN = './BSfiles/coilGxBS_Example_200626.biot';

[cont_xyz_Gx] = genisocont_sphcord( ver_mesh, fac_mesh, psi_mesh, xyz_ori, tz_bas, Niso,'BSFN',BSFN );

cont_Gx_all = cont_xyz_Gx;
save('./optresults/cont_Gx_all_Example_200626.mat','cont_Gx_all');

%% combine RO/shim coil with magnet model in single BS-file

% BSFN_mag = '/home/patmcd/Documents/HelmetMagnet/Magnet_v2/190327_proc/BSfiles/Magnet_190327_v01.biot';
% 
% BSFN_mag_ShimRO = 'BSfiles/Mag_ShimROy_190402_v04_test.biot';
% 
% BSfile_combine( BSFN_mag, BSFN, BSFN_mag_ShimRO );
