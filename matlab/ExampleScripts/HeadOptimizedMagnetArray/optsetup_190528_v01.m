%% filenames...
clear; close all;

fn_ROImask = '../../datafiles/ROImask_180809.mat';

%%% used in proc_ROImask_*** script
fn_MagROIsave = './precomp/MagROI_190528_v01.mat';

str_Dmat = '_190528_v01.mat';

fn_D1save = [ './precomp/D1' str_Dmat ];
fn_D3save = [ './precomp/D3' str_Dmat ];
fn_D5save = [ './precomp/D5' str_Dmat ];

fn_D1volsave = [ './precomp/D1vol' str_Dmat ];
fn_D3volsave = [ './precomp/D3vol' str_Dmat ];
fn_D5volsave = [ './precomp/D5vol' str_Dmat ];

fn_A_comp_save = [ './precomp/A_comp' str_Dmat ];

% fn_fnames = ['./precomp/fn_all' str_Dmat ];

%% generate ROI info
disp('generating ROI info...');
% proc_ROImask( fn_ROImask, fn_MagROIsave );
proc_ROImask_190528_v01( fn_ROImask, fn_MagROIsave );

%% generate D-matrices

disp('generating D1 matrix...');
proc_Dmatrix_L1( fn_MagROIsave, fn_D1save );
disp('generating D3 matrix...');
proc_Dmatrix_L3( fn_MagROIsave, fn_D3save );
disp('generating D5 matrix...');
proc_Dmatrix_L5( fn_MagROIsave, fn_D5save );

disp('generating D1vol matrix...');
proc_Dmatrix_vol_L1( fn_MagROIsave, fn_D1volsave );
disp('generating D3vol matrix...');
proc_Dmatrix_vol_L3( fn_MagROIsave, fn_D3volsave );
disp('generating D5vol matrix...');
proc_Dmatrix_vol_L5( fn_MagROIsave, fn_D5volsave );

%% generate A_comp matrix

disp('generating A_comp matrix...');
proc_A_comp( fn_MagROIsave, fn_A_comp_save );


%% generate J3, J5 (ie SH index permutation matrices)
gen_rotmat35_190214;