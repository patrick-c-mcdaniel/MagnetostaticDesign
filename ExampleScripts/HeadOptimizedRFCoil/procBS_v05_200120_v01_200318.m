fn_y_v05 = './BSfiles/BS_RF_y_zero_r08p5_200119_v05_dat_vol7mm.biot';

[ xx yy zz Bx By Bz Bm ] = BSdatread( fn_y_v05 );

SZdat = [61 43 43];

Bm_mat = reshape( Bm, SZdat );

xx_mat = reshape( xx, SZdat );
yy_mat = reshape( yy, SZdat );
zz_mat = reshape( zz, SZdat );

Bx_mat = reshape( Bx, SZdat );
By_mat = reshape( By, SZdat );
Bz_mat = reshape( Bz, SZdat );
Bxy_mat = sqrt( Bx_mat.^2 + By_mat.^2 );

mosaic_bg(permute(By_mat, [1 3 2]), 6, 7, 11, 'RF By; BS-simulation', [-1e-4 1e-4], [0 0 0], zeros(size(xx_mat)) ); colormap jet;

mosaic_bg(permute(Bx_mat, [1 3 2]), 6, 7, 12, 'RF Bx; BS-simulation', [-1e-4 1e-4], [0 0 0], zeros(size(xx_mat)) ); colormap jet;
mosaic_bg(permute(Bxy_mat, [1 3 2]), 6, 7, 21, 'RF Bxy; BS-simulation', [-1e-4 1e-4], [0 0 0], zeros(size(xx_mat)) ); colormap jet;
% mosaic_bg(permute(yy_mat, [1 3 2]), 6, 7, 22, 'y-coord; BS-simulation', [-0.147 0.147], [0 0 0], zeros(size(xx_mat)) ); colormap jet;
% mosaic_bg(permute(zz_mat, [1 3 2]), 6, 7, 23, 'z-coord; BS-simulation', [-0.147 0.147], [0 0 0], zeros(size(xx_mat)) ); colormap jet;
% 
% figure(11); mosaic_bg(permute(Bz_mat, [1 3 2]), 6, 7, 11, 'Shim Bz; BS-simulation', [-1e-3 3e-3], [0 0 0], zeros(size(xx_mat)) ); colormap jet;
% 
% figure(11); mosaic_bg(permute(Bz_mat, [1 3 2]), 6, 7, 11, 'Shim Bz; BS-simulation', [-1e-3 3e-3], [0 0 0], zeros(size(xx_mat)) ); colormap jet;

%% load ROI mask, etc.

dir_data = '/home/patmcd/Documents/HelmetMagnet/FieldMap/HelmetMagFM_200109/Data_200109_v01/';
dir_scr  = '/home/patmcd/Documents/HelmetMagnet/FieldMap/HelmetMagFM_200109/';
dir_precomp = '/home/patmcd/Documents/HelmetMagnet/Magnet_v2/190528_v01/precomp/';

% tmp = load('/home/patmcd/Documents/HelmetMagnet/Shim_Permanent/200110_meas/precomp/ShimROI_200110_v01.mat');

Mpc = load([ dir_precomp 'MagROI_190528_v01.mat' ]);
roi_ext = zeros( 18, 43, 43);
m_roi_ext = cat(1, Mpc.ROImsk_or_vol, roi_ext>0);

dispimgs( double(m_roi_ext), 'gray', 'axial');
% Bxyz_meas_ext = load('/home/patmcd/Documents/HelmetMagnet/FieldMap/proc_200109/proc_results/Bxyz_ROI_200110_v01.mat');

Bxy_ROI = Bxy_mat( m_roi_ext );

Bxy_mean = mean(Bxy_ROI(:));
figure(31); hold off;
histogram( Bxy_ROI );

%% ROI down towards feet

m_roi_feet1 = circshift( m_roi_ext, [2 0 0]);
m_roi_feet2 = circshift( m_roi_ext, [4 0 0]);

% m_roi_feet = (m_roi_feet1 | m_roi_feet2 ) & (~m_roi_ext);
m_roi_feet = xx_mat > 0.1 & (yy_mat.^2 + zz_mat.^2 <= 0.1^2 );

Bxy_zero = Bxy_mat( m_roi_feet );

Bxy_zero_mean = mean(Bxy_zero(:));
figure(31); hold on;
histogram( Bxy_zero, 'FaceColor', [1 0.2 0.2],'BinWidth',1e-6 );
grid on; grid minor;

xlim([0 1e-4]);
set(gcf,'color','white');

%% presentation figures 20/01/28
x1d = -0.147:0.007:0.273;
y1d = -0.147:0.007:0.147;
z1d = -0.147:0.007:0.147;


load('cm_red');

[ x3d y3d z3d ] = meshgrid( x1d, y1d, z1d);

figure(301); 
imagesc( z1d, x1d, squeeze( Bxy_mat(:,22,:) ) ); colormap jet; caxis([0 8e-5]); axis image;

figure(302); 
imagesc( z1d, x1d, squeeze( Bxy_mat(:,:,22) ) ); colormap jet; caxis([0 8e-5]); axis image;

load('/home/patmcd/Documents/HelmetMagnet/RF_B1LR/HelmetMagRF_200125/m_all_rot_ext_200128_v01.mat');

figure(311); set(gca,'Color','k');
imagesc( z1d, x1d, squeeze( Bxy_mat(:,22,:) ),'AlphaData', 0.3+0.7*squeeze(double(m_all_rot_ext(:,22,:)) )); colormap jet; caxis([0 8e-5]); axis image;
set(gcf,'color','white');


figure(312); set(gca,'Color',[0 0 0]);
imagesc( z1d, x1d, squeeze( Bxy_mat(:,:,22) ),'AlphaData', 0.3+0.7*squeeze(double(m_all_rot_ext(:,:,22)) )); colormap jet; caxis([0 8e-5]); axis image;
set(gcf,'color','white');
% 
% figure(311); 
% slice(x3d, y3d, z3d, Bxy_mat, [],[0],[0])

% mosaic_bg(permute(Bxy_mat, [1 3 2]), 6, 7, 21, 'RF Bxy; BS-simulation', [-1e-4 1e-4], [0 0 0], zeros(size(xx_mat)) ); colormap jet;