
%% load struct file

load('./optresults/MagStruct_OptExample_200625.mat','MagBlkSpec'    );
tht_adj     = MagBlkSpec.tht_adj;
phi_adj     = MagBlkSpec.phi_adj;
psi_adj     = MagBlkSpec.psi_adj;
xx_a        = MagBlkSpec.xx_a;
yy_a        = MagBlkSpec.yy_a;
zz_a        = MagBlkSpec.zz_a;
xsz_all     = MagBlkSpec.xsz_all;
ysz_all     = MagBlkSpec.ysz_all;
zsz_all     = MagBlkSpec.zsz_all;
Mm          = MagBlkSpec.Mm;

close all;


FACNUMs = ones(numel(xx_a),1);
FACNUMs(201:239) = 4;
FACNUMs(240:278) = 2;
FACNUMs(279:296) = 3;

FACNUMs(29:49)   = 4;
FACNUMs(50:70)   = 2;
FACNUMs(71:119)  = 3;
FACNUMs(120:140) = 2;
FACNUMs(141:161) = 4;

FACNUMs([ 26 166 ]) = 4;
% 
% FACNUMs([ 297:311 ]) = [ 3 3 5  2 2 2   0 0 0    4 4 2   1 4 5 ];
% FACNUMs([ 312:326 ]) = [ 3 3 6  2 2 2   0 0 0    4 4 2   1 4 6 ];
FACNUMs([ 297:311 ]) = [ 3 6 5  2 6 5   6 4 4    4 5 5   1 4 5 ];
FACNUMs([ 312:326 ]) = [ 3 5 6  2 5 6   5 4 4    4 6 6   1 4 6 ];

%% generate blocks (no tolerance; just to look at magnet array as STL)

TOLblk = 0.0e-3; % tolerance/design spacing between magnet blocks and coffers
HPLUG  = 0.00;   % Need at least depth of block (1" = 25.4mm)
% 
% [ ver_mag, fac_mag ] = Mspec2blocks_plug_booster( MagBlkSpec, TOLblk, 0, FACNUMs, 13 );
% [ ver_mag, fac_mag ] = Mspec2blocks_plug( MagBlkSpec, TOLblk, 0, FACNUMs, 12 );
[ ver_mag, fac_mag ] = Mspec2blocks_plug( MagBlkSpec, TOLblk, HPLUG, FACNUMs );
% 

fn_STL = './STLfiles/MagSTL_Blocks_200625.stl';
stlwrite( fn_STL, fac_mag, ver_mag );   

%% generate blocks (with extra spacing for block/coffer tolerance)

TOLblk = 0.2e-3; % tolerance/design spacing between magnet blocks and coffers
HPLUG  = 0.026;   % Need at least depth of block (1" = 25.4mm)
% 
% [ ver_mag, fac_mag ] = Mspec2blocks_plug_booster( MagBlkSpec, TOLblk, 0, FACNUMs, 13 );
% [ ver_mag, fac_mag ] = Mspec2blocks_plug( MagBlkSpec, TOLblk, 0, FACNUMs, 12 );
[ ver_mag, fac_mag ] = Mspec2blocks_plug( MagBlkSpec, TOLblk, HPLUG, FACNUMs );
% 

fn_STL = './STLfiles/MagSTL_BlockPlugs_200625.stl';
stlwrite( fn_STL, fac_mag, ver_mag );   


%% generate blocks (bounding blocks for coffers)

t_coffer = 3e-3;
EXTTOL = 0.5e-3;

[ ver_cof, fac_cof ] = Mspec2blocks_cof_exttol( MagBlkSpec, TOLblk+t_coffer, FACNUMs, EXTTOL );

fn_STL = './STLfiles/MagSTL_Coffers_200625.stl';
stlwrite( fn_STL, fac_cof, ver_cof );   

%  view([0 90 0]);