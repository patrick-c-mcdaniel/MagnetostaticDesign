close all;
clear;

% load('/home/patmcd/Documents/HelmetMagnet/Magnet_v2/180926_v01/Mout_20180926T193744.mat','Mout');
% dir_rem = '/home/patmcd/Documents/HelmetMagnet/Magnet_v2/190523_v01/';
%% run preparation scripts
% ShimPerm_ptsROI_190626_v01;
% ShimPerm_genD_190626_v01;

%% load D-matrices for optimization field computation (on surface)
% load(fn_D1save,'D2d');
% load(fn_D3save,'D3_2d');
% load(fn_D5save,'D5_2d');

load('./precomp/D1comp_200110_v01.mat');
load('./precomp/ShimROI_200110_v01.mat');

% %% load J-matrices (coordinate permutation matrices in SH-bases)
% load([  '../SharedMat/J3.mat']);
% load([  '../SharedMat/J5.mat']);

% %% load M-symmetry matrix A_comp
% load(fn_A_comp_save);

% %% load old solution and target ROI coordinates, mask, etc.
% % load('/home/patmcd/Documents/HelmetMagnet/Magnet_v2/180926_v01/Mout_20180926T193744.mat','Mout');
% load(fn_MagROIsave);


%%% initial guess



%% initialize fmincon

%% Optimization parameters; contraint values; variable options (discrete set for ga)
BMEAN = 0.07;

Mred = zeros(Nmag, 1);

lbM = 0*ones(size(Mred));   % Optimization variables can take the values {0, 1, 2}. During field computation, each parameter is multiplied 
ubM = 2*ones(size(Mred));   %   by Mdel (the magnet block "discretization step" size; has units of magnetication) to give a physically-meaningful
                            %   quantity.
MUz = 4*pi*1e-7;
% Br = 1.32;
Br = 1.38; % N45 material

Mdel =  2*1/MUz*Br*(1/8)^3*(0.0254^3);


%% Define functions for fmincon (cost, constraint, plot)

%%% cost function
% fcost = @(Mred) ( cost_Buniform_nomean( Bcomp135_sz3_tps_1var( D135, J3, J5, Mxyz_to_sz3tps(Mred2Mxyz(Mred, inds_mid, inds_pos, inds_neg),tps(:,3)), errco_all ) ) );
fcost = @(Mred) ( cost_Buniform_nomean_Btarg( Mdel*D1comp*Mred'+Btarg ) );
% fcost = @(Mred) ( cost_BzL2_Bxy0( Bcomp135_Mred( D135red, J3, J5, Mred, sssred, Nmid, errco_all ) ) );


%%% nlcon function
% fnlcon3 = @(Mred, Mmax, BMEAN) deal( cat(1, 5e4*nlcon_Bmin( Bcomp135_sz3_tps_1var( D135, J3, J5, ...
%                                                                 Mxyz_to_sz3tps(Mred2Mxyz(Mred, inds_mid, inds_pos, inds_neg),tps(:,3)), ...
%                                                                 errco_all ), ...
%                                                             BMEAN), ...
%                                                 nlcon_Mmax(Mred2Mxyz(Mred, inds_mid, inds_pos, inds_neg), Mmax)), [ ] ); 
fnlcon3 = @(Mred, BMEAN) ( nlcon_Bmin_Btarg( Mdel*D1comp*Mred'+Btarg, BMEAN)); 

fnlcon  = @(Mred) deal( [fnlcon3(Mred, BMEAN)], []);

%%% plot function (each iteration)
load('cm_rb.mat');

% for plot, use the full 3D volume field computation
fplot   = @(optinfo, state, flag) ga_plot_helper(optinfo, state, flag, Mdel, D1comp_vol, Btarg_vol, ROImsk_or_vol, cm_rb);
% fplot   = @(optinfo, state, flag) test_old_ga_plot(optinfo, state, flag);


% Mxyz = sz3tps_to_Mxyz( [sz3 tps] );
% Mxyz(1:20:end)=0;

%%% set fmincon options
%
% Note: names for optfmincon fields+conventions for values are different
% for R2013a (wtf...)

% optga = optimoptions('ga');
% optga = gaoptimset();
% optga.MaxGenerations             = 250000;
% optga.MaxStallGenerations        = 100;
% optga.Display            = 'iter';
% % optga.OutputFcn          = fplot;
% % optga.PlotFcn            = fplot;
% optga.PlotFcns            = fplot;
% optga.PlotInterval       = 20;

%% optimization setup parameters; might need to tweak for R2015a/R2013a/R2017a Matlab (sigh...)

optga = gaoptimset('PopulationSize', 150, 'Display', 'diagnose' ,'PlotFcns',{ fplot },'StallGenLimit',100,...
    'Generations', 2500, 'PlotInterval', 20, 'TolFun', 1e-10);  
optga_save = optga;

%%% the plot function includes fplot, which uses the D-matrices in its definition, which in turn are very large
%%% removing the plotfcns thus lets you save the optimization settings without taking up tons of space
optga_save.PlotFcns = {};

IntConMask = 1:Nmag;

AA = [];
bb = [];

%% run ga; save the result from each iteration in ./optresults

% I normally would run ~25 iterations; I set this to 3 just because it'll run more quickly
Nrun = 3;

for irr = 1:Nrun
    tic
    [ Mout ] = ga( fcost, Nmag, [], [], [], [], lbM', ubM', fnlcon, IntConMask, optga );
    toc
    T_run = toc;

% 
% tic
% [ Mout ] = fmincon( fcost, Mred, [], [], [], [], lbM, ubM, fnlcon, optfmincon );
% toc
% T_run = toc;


% strout = datestr(datetime,30);
% 
% save( ['Mout_' strout '.mat'], 'Mout','T_run','optfmincon','fcost','Mmax','Bmin' );



    strout = datestr(now,30);
    MachineName = 'jotunn.mit.edu';
    save( ['./optresults/Mout_' num2str(irr,'%03u') '_' strout '.mat'], 'Mout','T_run','optga_save','BMEAN','MachineName');
end
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
