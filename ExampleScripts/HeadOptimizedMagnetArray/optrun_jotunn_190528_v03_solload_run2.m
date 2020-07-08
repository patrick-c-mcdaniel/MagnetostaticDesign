% close all;
% clear;

% load('/home/patmcd/Documents/HelmetMagnet/Magnet_v2/180926_v01/Mout_20180926T193744.mat','Mout');
% dir_rem = '/home/patmcd/Documents/HelmetMagnet/Magnet_v2/190523_v01/';

%% load D-matrices for optimization field computation (on surface)
load(fn_D1save,'D2d');
load(fn_D3save,'D3_2d');
load(fn_D5save,'D5_2d');

%% load J-matrices (coordinate permutation matrices in SH-bases)
load([  './precomp/J3.mat']);
load([  './precomp/J5.mat']);

%% load M-symmetry matrix A_comp (used to reduce number of computations by exploiting magnet symmetry)
load(fn_A_comp_save);

%% load target ROI coordinates, mask, initial magnet design, etc. 
% load('/home/patmcd/Documents/HelmetMagnet/Magnet_v2/180926_v01/Mout_20180926T193744.mat','Mout');

load(fn_MagROIsave);


%%% initial guess

Mxyz = mag2vec( mx_ba, my_ba, mz_ba )/2;

% 
% Nmag = numel(xx_a);
% 
% % Mxyz_full = zeros(Nmag*3,1);
% % 
% Mxyz_296_obj = load('/home/patmcd/Documents/HelmetMagnet/Magnet_v2/190523_v01/Mout_tmp_save_190524_122141.mat','Mout');
% 
% Nlg  = numel(Mxyz_296_obj.Mout);
% 
% Mxyz_326_obj = load('./Mout_20190529T004557.mat','Mout');
% Mred_326 = Mxyz_326_obj.Mout;
% 
% Mxyz_326 = Mred2Mxyz( Mred_326, inds_mid, inds_pos, inds_neg );
% % 
% Mxyz(1:Nmag*3) = Mxyz_326(1:Nmag*3);
% Mxyz = load('./optresults/Mout_tmp_save_190530_064431.mat','Mout');

%%% The number of total magnet blocks + the number of "large" magnet blocks (used for setting block size constraints)
Nmag = 326; 
Nlg  = 296;
% Mxyz = Mxyz_full;

%%%% initial go (use analytically-generated "halbach" design")
Mred = Mxyz2Mred( Mxyz, inds_mid, inds_pos );

%%% run2 initialization (use magnet design solution from earlier iteration of optimization; increase block sizes by 10%) 
Mout_run1 = load('./optresults/Mout_20200623T041906_Run1Example.mat');
Mred = Mout_run1.Mout * 1.1;

%% load matrices for 3D volume field computation (for visualization)
D1vol = load(fn_D1volsave,'D2d');
D3vol = load(fn_D3volsave,'D3_2d');
D5vol = load(fn_D5volsave,'D5_2d');

D135vol = [D1vol.D2d, D3vol.D3_2d, D5vol.D5_2d];
clear D1vol D3vol D5vol;

%% load error coefficients for block SH terms
dn_err = '../BlockMultipoleComp/matfiles/';
fn_err_L1M0 = 'pf_1p05_L1M0.mat'; load([dn_err fn_err_L1M0]);
fn_err_L3M0 = 'pf_1p05_L3M0.mat'; load([dn_err fn_err_L3M0]);
fn_err_L3M2 = 'pf_1p05_L3M2.mat'; load([dn_err fn_err_L3M2]); 
fn_err_L5M0 = 'pf_1p05_L5M0.mat'; load([dn_err fn_err_L5M0]);
fn_err_L5M2 = 'pf_1p05_L5M2.mat'; load([dn_err fn_err_L5M2]);
fn_err_L5M4 = 'pf_1p05_L5M4.mat'; load([dn_err fn_err_L5M4]);

errco_all = cat(1, pf_1p05_L1M0, ...
                   pf_1p05_L3M0, ...
                   pf_1p05_L3M2, ...
                   pf_1p05_L5M0, ...
                   pf_1p05_L5M2, ...
                   pf_1p05_L5M4);
               
%% test using final comsol model from 181104 (incl. psi angles)
% load('../../datafiles/MagSpec_181104_v01.mat','tht_adj','phi_adj','psi_adj','Mm','xsz_all','ysz_all','zsz_all'    );
% 
% phi_adj([1:7 92:98]) = 90*[ -1 1 1 1 -1 -1 1 1 -1 -1 1 1 -1 -1];
% 
% sz3 = [xsz_all(:) ysz_all(:) zsz_all(:)];
% tps = [tht_adj(:), phi_adj(:), psi_adj(:)]*pi/180;
% 
% 
% tpsfull = zeros( Nmag, 3);
% tpsfull( 1:size(tps,1), :) = tps;
% tpsfull( (Nmag-29):(Nmag-15), 3 ) = +pi/2;
% tpsfull( (Nmag-14):(Nmag+0 ), 3 ) = -pi/2;
% 
% tps = tpsfull;


%% configure compressed matrices; etc

%%% D-matrix
D135 = cat(2, D2d, D3_2d, D5_2d);

D135red = D135 * A_comp;

%% set block "psi" angles to zero (not optimized over)

sss = zeros(Nmag, 1); %tps(:,3);
sssred = cat(1, sss( inds_mid(:)), sss(inds_pos(:)) );

Nmid = numel(inds_mid);

%% set up Btarg

% shimmat = load('/home/patmcd/Documents/HelmetMagnet/Magnet_v2/ShimCoil/Shim_Comsol181008/ShimOptSetup_181008_v01.mat');
% 
% Btarg_mat = zeros(size(shimmat.ROImsk_or));
% Btarg_mat( shimmat.ROImsk_or ) = shimmat.Btarg;
% Btarg_mat_zfl = flip(Btarg_mat,3);
% Btarg = Btarg_mat_zfl( flip(shimmat.ROImsk_or,3) );
% Btarg = cat(1, zeros(size(Btarg')), zeros(size(Btarg')), Btarg');
% Btarg = Btarg(:);

%% initialize fmincon


%%%%%%%%%%%%% Configure block size constraints

MUZ  = 4*pi*1e-7;               % SI units

%%% Mmax/Bmean for constraints
Mmax_ID_1 = 1.42 * (0.0254)^3 * 1*1*1 / MUZ;    % SI units
Mmax_ID_2 = 1.42 * (0.0254)^3 * 0.5*1*1 / MUZ;    % SI units

Mmax = Mmax_ID_1;

BMEAN = 0.07;                    % unit [T]

%%% upper and lower bounds on Mx, My, Mz
lbM = -1*Mmax_ID_1*ones(numel(mm_ba), 3);
lbM(end-29:end) = -1*Mmax_ID_2;
ubM =  1*Mmax_ID_1*ones(numel(mm_ba), 3);
ubM(end-29:end) =  1*Mmax_ID_2;

lbM = lbM(:);
ubM = ubM(:);
% 
% Mmax_all = zeros(Nmag*3,1);
% Mmax_all(1:3*Nlg) = Mmax_ID_1;
% Mmax_all((3*Nlg+1):end) = Mmax_ID_2;

Mmax_all = zeros(Nmag,1);
Mmax_all(1:Nlg) = Mmax_ID_1;
Mmax_all((Nlg+1):end) = Mmax_ID_2;

%%%%%%%%%%%% Define functions needed for fmincon

%%% cost function
% fcost = @(Mred) ( cost_Buniform_nomean( Bcomp135_sz3_tps_1var( D135, J3, J5, Mxyz_to_sz3tps(Mred2Mxyz(Mred, inds_mid, inds_pos, inds_neg),tps(:,3)), errco_all ) ) );
fcost = @(Mred) ( cost_Buniform_nomean( Bcomp135_Mred( D135red, J3, J5, Mred, sssred, Nmid, errco_all ) ) );
% fcost = @(Mred) ( cost_BzL2_Bxy0( Bcomp135_Mred( D135red, J3, J5, Mred, sssred, Nmid, errco_all ) ) );


%%% nlcon function
% fnlcon3 = @(Mred, Mmax, BMEAN) deal( cat(1, 5e4*nlcon_Bmin( Bcomp135_sz3_tps_1var( D135, J3, J5, ...
%                                                                 Mxyz_to_sz3tps(Mred2Mxyz(Mred, inds_mid, inds_pos, inds_neg),tps(:,3)), ...
%                                                                 errco_all ), ...
%                                                             BMEAN), ...
%                                                 nlcon_Mmax(Mred2Mxyz(Mred, inds_mid, inds_pos, inds_neg), Mmax)), [ ] ); 
fnlcon3 = @(Mred, Mmax, BMEAN) deal( cat(1, 5e4*nlcon_Bmin( Bcomp135_Mred( D135red, J3, J5, Mred, sssred, Nmid, errco_all ), BMEAN), ...
                                                nlcon_Mmax_vec(Mred2Mxyz(Mred, inds_mid, inds_pos, inds_neg), Mmax_all)), [ ] ); 

fnlcon  = @(Mred) fnlcon3(Mred, Mmax_ID_1, BMEAN);

%%% plot function (each iteration)
%%%   The function also saves current fmincon state (including solution at iteration) in case it crashes or you manually want to end the solver
%%%   prematurely.
load('cm_rb.mat');

% for plot, use the full 3D volume field computation
fplot   = @(Mred, optinfo, state) plot_mag_bfield_L135(Mred2Mxyz(Mred, inds_mid, inds_pos, inds_neg), optinfo, state, ...
                                                       Bcomp135_sz3_tps_1var( D135vol, J3, J5, ...
                                                            Mxyz_to_sz3tps(Mred2Mxyz(Mred, inds_mid, inds_pos, inds_neg),sss), ...
                                                            errco_all ), ...
                                                       ROImsk_or_vol, xx_a, yy_a, zz_a, cm_rb, 1 );



% Mxyz = sz3tps_to_Mxyz( [sz3 tps] );
% Mxyz(1:20:end)=0;

%%%%%%%%%%%%%% set fmincon options
%
% Note: names for optfmincon fields+conventions for values are different
% for R2013a, R2015a (wtf...)

%%% R2017a
% optfmincon = optimoptions('fmincon');
% optfmincon.MaxIter            = 100000;
% optfmincon.MaxFunEvals        = 10000000;
% optfmincon.UseParallel        = 0; %1;
% optfmincon.Display            = 'iter';
% optfmincon.Algorithm          = 'interior-point';
% optfmincon.OutputFcn          = fplot;
% optfmincon.StepTolerance      = eps;

%%% R2015a
optfmincon = optimoptions('fmincon');
optfmincon.MaxIter            = 100000;
optfmincon.MaxFunEvals        = 10000000;
optfmincon.UseParallel        = 0; %1;
optfmincon.Display            = 'iter';
optfmincon.Algorithm          = 'interior-point';
optfmincon.OutputFcn          = fplot;
optfmincon.TolX               = eps;

%%% the plot function includes fplot, which uses the D-matrices in its definition, which in turn are very large
%%% removing the plotfcns thus lets you save the optimization settings without taking up tons of space
optfmincon_save = optfmincon;
optfmincon_save.OutputFcn = {};

% return
%% run fmincon

tic
[ Mout ] = fmincon( fcost, Mred, [], [], [], [], lbM, ubM, fnlcon, optfmincon );
toc
T_run = toc;

% strout = datestr(datetime,30);
% 
% save( ['Mout_' strout '.mat'], 'Mout','T_run','optfmincon','fcost','Mmax','Bmin' );



strout = datestr(now,30);
MachineName = 'jotunn.mit.edu';
save( ['optresults/Mout_' strout '.mat'], 'Mout','T_run','optfmincon_save','Mmax','BMEAN','MachineName');
%% plot resulting field map

% Bout = D2d * Mout;
% 
% [ Bx By Bz ] = vec2field( Bout );
% 
% Bm = sqrt( Bx.^2 + By.^2 + Bz.^2 );
% 
% Bm_mat = zeros(size(ROImsk_or));
% Bm_mat(ROImsk_or) = Bm;
% 
% mosaic_bg( permute(Bm_mat, [1 2 3]), 5, 9, 11, [], [], [0 0 0],permute(ROImsk_or, [1 2 3])==0 ); colormap jet; colorbar; caxis([min(Bm(:)) max(Bm(:))]);
% mosaic_bg( permute(Bm_mat, [2 3 1]), 5, 9, 12, [], [], [0 0 0],permute(ROImsk_or, [2 3 1])==0 ); colormap jet; colorbar; caxis([min(Bm(:)) max(Bm(:))]);
% mosaic_bg( permute(Bm_mat, [3 1 2]), 5, 9, 13, [], [], [0 0 0],permute(ROImsk_or, [3 1 2])==0 ); colormap jet; colorbar; caxis([min(Bm(:)) max(Bm(:))]);
