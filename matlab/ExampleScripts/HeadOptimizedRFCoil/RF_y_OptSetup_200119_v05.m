%% Notes 20/01/19
%    - y-polarized field
%    - R_roi = 8.5cm (coil is narrow)


%% needed for coil optimization:
%       - set of target points
%       - field (Bx, By, Bz) at all target points
%       - mesh for designing stream function
%       - regularization term
%       - regularization preconditioning (weighting) matrix


%% load files as needed

%%% target point coordinates; in:  x3tar, y3tar, z3tar
% surface mask in: ROImsk_or
% volume mask in:  ROImsk_or_vol
load('../HeadOptimizedMagnetArray/precomp/MagROI_190528_v01.mat');

% fn_coil_red = 'coil_red_RF_y_zero_r08p5_200119_v05.mat';
%%% B-field; load BS-simulation file 
% dn_bsdat = '/home/patmcd/Documents/HelmetMagnet/Magnet_v2/190327_proc/BSfiles/';
% fn_bsdat = 'Magnet_190327_v01_datXYZ.biot';
% [ xxBs, yyBs, zzBs, BxBs, ByBs, BzBs, BmBs ] = BSdatread([ dn_bsdat fn_bsdat ]);

%%% stream function/coil surface mesh
[listNode, listTriangle, tri] = importMeshWavefront('./mesh/mesh_interp_200119_v01.obj');


coil.listNode = listNode;
coil.listTriangle = listTriangle;
coil.tri = tri;

% clear coil;
% load coil_190402_v02.mat;

%% generate target field/ROI

R_ROI = 0.08;

SZ3d = [43 43 43];

xx_mat = x3tar; %zeros(size(ROImsk_or_vol));
yy_mat = y3tar; 
zz_mat = z3tar; 
r3tar = sqrt( x3tar.^2 + y3tar.^2 + z3tar.^2 );
ryz3tar = sqrt( y3tar.^2 + z3tar.^2 );
% Bm_mat_all = reshape(BmBs, SZ3d);

x_ROI = x3tar(ROImsk_or_vol);
xmax = max(x_ROI(:));

% ROImsk_or_vol_xext = ROImsk_or_vol | circshift(ROImsk_or_vol,[3 0 0]) | circshift(ROImsk_or_vol,[6 0 0]);
ROImsk_targ_unity = (x3tar<=xmax) & ((x3tar>=0 & ryz3tar<R_ROI) | r3tar<R_ROI);
ROImsk_targ_zeros = (x3tar>xmax) & ((x3tar>=0 & ryz3tar<R_ROI) | r3tar<R_ROI);

% By_tar_mat = zeros(size(xx_mat));
% By_tar_mat(ROImsk_targ) = 1;

%%% unity points
x1t = xx_mat(ROImsk_targ_unity);
y1t = yy_mat(ROImsk_targ_unity);
z1t = zz_mat(ROImsk_targ_unity);

By1t = ones(size(x1t));

%%% zero points
x0t = xx_mat(ROImsk_targ_zeros);
y0t = yy_mat(ROImsk_targ_zeros);
z0t = zz_mat(ROImsk_targ_zeros);

By0t = zeros(size(x0t));

xt_vec = xx_mat(:);
yt_vec = yy_mat(:);
zt_vec = zz_mat(:);

ROImsk_tar = ROImsk_or_vol;


pmat = [3 1 2];

mosaic(permute( double(ROImsk_targ_unity), pmat) , 5, 7, 106); colormap jet; colorbar; caxis([-5e-4 5e-4]); title('unity points');
mosaic(permute( double(ROImsk_targ_zeros), pmat) , 5, 7, 107); colormap jet; colorbar; caxis([-5e-4 5e-4]); title('zero points');

%%% x,y,z target ROI coords; Btarg

Btarg = cat(1, By1t(:), By0t(:));
xROI = cat(1, x1t(:), x0t(:));
yROI = cat(1, y1t(:), y0t(:));
zROI = cat(1, z1t(:), z0t(:));

%% plot of mesh with target ROI points

figure(111); hold off;
patchsexy( 'Faces', listTriangle, 'Vertices', listNode, 'FaceAlpha', 0.6 ); hold on;

plot3( x1t, y1t, z1t, 'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.8 0 0] );
plot3( x0t, y0t, z0t, 'o', 'MarkerSize', 5, 'MarkerFaceColor', [0 0.8 0 ]  );

%% create target field for shim optimization

% %%% use z-component of B0 as target field, by convention
% Btarg = Bt_mat(ROImsk_or_vol);
% Btarg_mat = Bt_mat;
% 
% % mosaic(permute( Btarg_mat, pmat) , 5, 7, 111); colormap jet; colorbar; caxis([min(Btarg(:)) max(Btarg(:))]); title('Shim coil target B-field [T]');
% 
% %%% add Gy gradient to target field
% 
% Gro_all = 0; %(0:4:20)*1e-3;    % units: T/m
% 
% 
% Ngro = numel(Gro_all);
% 
% pYZ = [2 3 1];
% PXZ = [3 1 2];
% 
% 
% for igr = 1:Ngro
% 
%     Bry = yt_vec * Gro_all(igr);
%     
%     Btarg = Btarg + Bry;
% end
% 
% save('ShimOptSetup_190402_v04.mat','Btarg','ROImsk_tar','x3tar','y3tar','z3tar');

%% Import meshing
% If reduction = 1, the matrix system will be reduced, in order to set all
% the border node to the same value, in order to respect the divergence free criteria
% i.e. the node vector will be reorgenized in order to have all the border
% node on top, in order to facilitate the reduction of all the system

% [coil.listNode,coil.listTriangle,coil.tri] = importMeshWavefront('./mesh/ShimSurf_v02_180925_v06obj.obj');
coil.center = [0 0 0];
coil.reduction = 1;
coil.rateIncreasingWire = 1;

%% Target point definition

% degreeMax = 21;
% orderMax = 21;
% rhoReference = 0.05; % radius of reference
% rk1 = createTargetPointGaussLegendreAndRectangle_offset_SS(rhoReference,degreeMax,orderMax);
% 
% rk1 = rk1( (rk1(:,1)>=-0.04)&(rk1(:,1)<=0.04)& ...
%          ((rk1(:,2).^2 + rk1(:,3).^2)<=(0.05^2)) ,:);
%      
% rk2 = createTargetPointGaussLegendreAndRectangle_offset_SS(rhoReference*2/3,15,15);
% 
% rk2 = rk2( (rk2(:,1)>=-0.04)&(rk2(:,1)<=0.04)& ...
%          ((rk2(:,2).^2 + rk2(:,3).^2)<=(0.05^2)) ,:);
%      
% rk3 = createTargetPointGaussLegendreAndRectangle_offset_SS(rhoReference*1/3,9,9);
% 
% rk3 = rk3( (rk3(:,1)>=-0.04)&(rk3(:,1)<=0.04)& ...
%          ((rk3(:,2).^2 + rk3(:,3).^2)<=(0.05^2)) ,:);
%      
%      rk = cat(1, rk1, rk2, rk3);
% 
% xROI = x3tar(ROImsk_tar);
% yROI = y3tar(ROImsk_tar);
% zROI = z3tar(ROImsk_tar);

rk = cat(1, xROI(:)', yROI(:)', zROI(:)');
rk = rk';
%% Then we habe to calculate the field in a given direction

% Initialize the target ampltiude
% bc(1).coefficient = zeros(degreeMax+1,orderMax+1);
% bs(1).coefficient = zeros(degreeMax+1,orderMax+1);
% bc(2).coefficient = zeros(degreeMax+1,orderMax+1);
% bs(2).coefficient = zeros(degreeMax+1,orderMax+1);
% bc(3).coefficient = zeros(degreeMax+1,orderMax+1);
% bs(3).coefficient = zeros(degreeMax+1,orderMax+1);

% MRI Gradient dBz/dy
targetCoil = 'By';
% bs(3).coefficient(2,2) = 0.005*rhoReference; % Drive
% 
% B  = RebuildField7bis(bc,bs,rhoReference, rk,'sch');
% coil.btarget = [B(3,:)];
coil.btarget = Btarg(:);
%DisplayFieldOnSphere( B,rk,'TargetField' )
% 
% clear('B');
 
 %% In order to calculate the resistance matrix, we have to provide some data :

% For the coil
coil.wireThickness = 0.001; % (meter) Thickness of the conductor
coil.wireWidth = 0.002; % (meter) Thickness of the conductor
coil.wireSurface = coil.wireThickness*coil.wireWidth; % in meter %5mmx5mm is equivalent to the number used in Timo's coil or the 7.5*7.5 litz wire
coil.fillFactor = 1;
coil.rhoCopper = 1.68*10^-8; % (Ohm*m) resistivity of the copper
coil.wireResistivity = coil.rhoCopper/coil.fillFactor;  % (Ohm*m) resistivity of the wire

%% Part to calculate
%optimizationType = 'QP';
%coil.error = 0.05;
optimizationType = 'standardTikhonov';
% reg=4.0e-6; % initial value from code
reg=6.0*10^-6;
calculateR = 1;
calculateL = 1;
calculateLwp = 0;

coil.startingWireNumber = 2;
coil.distanceBetween2Wire = 20*10^-3;
coil.rateIncreasingWire = 2;

%% Set the variable for the Field calculation

% X position in meter
coil.x_Start =-0.120;%0.136
coil.x_Stop = 0.120;
coil.x_Step = 0.01;
coil.x_Value = coil.x_Start:coil.x_Step:coil.x_Stop;
% Y position in meter
coil.y_Start = -0.120;
coil.y_Stop = 0.120;
coil.y_Step = 0.01;
coil.y_Value = coil.y_Start:coil.y_Step:coil.y_Stop;
% Z position in meter
coil.z_Start =-0.002;
coil.z_Stop = 0.002;
coil.z_Step = 0.002;
coil.z_Value = coil.z_Start:coil.z_Step:coil.z_Stop;

coil.current = 1;
coil.sphere_radius = 0.025;
coil.coil_radius = 0.1;
coil.coil_length = 0.2;

save('./precomp/OptSetupRFLR_Example_200626.mat','Btarg','ROImsk_tar','x3tar','y3tar','z3tar','optimizationType','reg','calculateR','calculateL','calculateLwp','rk','targetCoil');

save('./precomp/coil_OptSetupRFLR_Example_200626.mat','coil');