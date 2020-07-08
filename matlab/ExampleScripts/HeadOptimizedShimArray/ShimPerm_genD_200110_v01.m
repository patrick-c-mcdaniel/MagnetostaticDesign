close all; clear;

load( './precomp/ShimROI_200110_v01.mat');

fn_MagROI = './precomp/ShimROI_200110_v01.mat';

%% create full D1 matrix ( (3*Nm)x(3*Nb) )

[ D1full ] = proc_Dmatrix_L1( fn_MagROI, [] );

%% create full volume D1 matrix ( (3*Nm)x(3*Nbvol) )

[ D1full_vol ] = proc_Dmatrix_vol_L1( fn_MagROI, [] );

%% create M-compression matrix

Rmx = diag(mx_ba);
Rmy = diag(my_ba);
Rmz = diag(mz_ba);
Rm_3d = cat(3, Rmx, Rmy, Rmz);

Rm_3d = permute( Rm_3d, [3 1 2] );
Rm_2d = reshape( Rm_3d, [ 3*Nmag, Nmag] );

figure(21); imagesc( Rm_2d ); axis image;

%% create B-compression matrix

%%%%%% edit once COMSOL simulations have run
% fn_Btarg = '/home/patmcd/Documents/HelmetMagnet/Magnet_v2/190528_v01/Comsol/Datafiles/Dat_190528_v03_discrete_N023_w2_Br1p42_mur1p05_20190702T105618_vol7m.txt';
% 
% dat_comsol = dlmread(fn_Btarg);
% 
%     SZCM = [43 43 43];
% 
%     Bx_cm = reshape(dat_comsol(:,4), SZCM);
%     By_cm = reshape(dat_comsol(:,5), SZCM);
%     Bz_cm = reshape(dat_comsol(:,6), SZCM);
%     Bm_cm = sqrt( Bx_cm.^2 + By_cm.^2 + Bz_cm.^2 );
%  
% Btarg = Bm_cm(ROImsk_or);
% Bxtarg = Bx_cm(ROImsk_or);
% Bytarg = By_cm(ROImsk_or);
% Bztarg = Bz_cm(ROImsk_or);

Bxyz_meas_ext = load('../../datafiles/Bxyz_ROI_200110_v01.mat');

Bx_cm = zeros(size(Bxyz_meas_ext.Bm_ROI_ext_mat));
By_cm = zeros(size(Bxyz_meas_ext.Bm_ROI_ext_mat));
Bz_cm = Bxyz_meas_ext.Bm_ROI_ext_mat ;
Bm_cm = sqrt( Bx_cm.^2 + By_cm.^2 + Bz_cm.^2 );
    
Btarg = Bm_cm(ROImsk_or);
Bxtarg = Bx_cm(ROImsk_or);
Bytarg = By_cm(ROImsk_or);
Bztarg = Bz_cm(ROImsk_or);

    

nBx = Bxtarg./Btarg;
nBy = Bytarg./Btarg;
nBz = Bztarg./Btarg;

Rbx = diag(nBx);
Rby = diag(nBy);
Rbz = diag(nBz);

Rb_3d = cat(3, Rbx, Rby, Rbz);
Rb_2d = reshape( permute( Rb_3d, [1 3 2] ), [Ntar, 3*Ntar] );

figure(31); imagesc( Rb_2d ); axis image;

%% create B-compression matrix (vol)
 
Btarg_vol = Bm_cm(ROImsk_or_vol);
Bxtarg_vol = Bx_cm(ROImsk_or_vol);
Bytarg_vol = By_cm(ROImsk_or_vol);
Bztarg_vol = Bz_cm(ROImsk_or_vol);

nBx_vol = Bxtarg_vol./Btarg_vol;
nBy_vol = Bytarg_vol./Btarg_vol;
nBz_vol = Bztarg_vol./Btarg_vol;

Rbx_vol = diag(nBx_vol);
Rby_vol = diag(nBy_vol);
Rbz_vol = diag(nBz_vol);

Rb_3d_vol = cat(3, Rbx_vol, Rby_vol, Rbz_vol);
Rb_2d_vol = reshape( permute( Rb_3d_vol, [1 3 2] ), [Ntar_vol, 3*Ntar_vol] );

% figure(31); imagesc( Rb_2d ); axis image;

%% compress D1-matrix

D1comp = Rb_2d * D1full * Rm_2d;
D1comp_vol = Rb_2d_vol*D1full_vol*Rm_2d;

fn_D1comp = './precomp/D1comp_200110_v01.mat';
save(fn_D1comp, 'D1comp','D1full','Rm_2d','Rb_2d','D1comp_vol','Bm_cm','Btarg','Btarg_vol');

