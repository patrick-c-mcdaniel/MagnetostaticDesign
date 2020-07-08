%% Filenames to process

fn_Mout           = './optresults/Mout_Run2Example.mat';       % Solution to optimization 
fn_MagROI         = './precomp/MagROI_190528_v01.mat';  % ROI/magnet geometry information

% Pre-existing magnet array spec; only used to obtain "psi"-angles for the
% magnet blocks (these angles are not determined by the optimization
% solution, and must be specified on their own). 
%
% Note that this file only needs one variable for the purposes of this
% script:
%   psi_adj     : Nblk x 1 array; values of block psi angles (in DEGREES)
%
fn_MagBlkSpec     = '../../datafiles/Mag_psi_adj_ExampleNonzero_200625.mat';   
fn_MagBlkSpec_mod = '../../datafiles/Mag_psi_adj_ExampleNonzero_200625_mod.mat';   

%% manually modify psi angles
% 
% load(fn_MagBlkSpec)
% 
% psi_adj(12) = -45; 
% psi_adj(13) = -50;
% psi_adj(14) = 0;
% psi_adj(20) = -20;
% psi_adj(21) = 0;
% 
% %%%%
% psi_adj(180) = 45;
% psi_adj(181) = 50;
% psi_adj(182) = 0;
% psi_adj(174) = 20;
% psi_adj(175) = 0;
% 
% %%%%
% psi_adj(297) = -45;
% psi_adj(298) = -30;
% psi_adj(299) = 5;
% 
% psi_adj(300:302) = -30;
% 
% psi_adj(306:308) = 0;
% 
% psi_adj(309) = 90;
% psi_adj(310) = -5;
% psi_adj(311) = 0;
% 
% psi_adj(312:314) = 0;
% 
% %%%%
% psi_adj(312) = 45;
% psi_adj(313) = 30;
% psi_adj(314) = -5;
% 
% psi_adj(315:317) = 30;
% 
% psi_adj(321:323) = 0;
% 
% psi_adj(324) = -90;
% psi_adj(325) = 5;
% psi_adj(326) = 0;
% 
% save(fn_MagBlkSpec_mod,'psi_adj');

%% Output filenames
% [ ver_all, fac_all ] = Mred2blocks( fn_Mout, fn_MagROI, fn_MagBlkSpec );

fn_BSsave = './BSfiles/BS_example_200625.biot';           % Output filename at which to save Biot-Savart file
fn_MagStruct = './optresults/MagStruct_OptExample_200625.mat';          % Output filename at which to save the magnet array structure


%% Function that saves a B-S file and returns a magnet array structure
[ MagBlkSpec ] = Mred2BSfile( fn_Mout, fn_MagROI, fn_MagBlkSpec_mod, fn_BSsave );


save( fn_MagStruct, 'MagBlkSpec' );