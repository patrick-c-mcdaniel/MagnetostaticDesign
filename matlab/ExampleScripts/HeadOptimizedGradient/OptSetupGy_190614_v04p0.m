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

% %%% B-field; load BS-simulation file 
% dn_bsdat = '/home/patmcd/Documents/HelmetMagnet/Magnet_v2/190327_proc/BSfiles/';
% fn_bsdat = 'Magnet_190327_v01_datXYZ.biot';
% [ xxBs, yyBs, zzBs, BxBs, ByBs, BzBs, BmBs ] = BSdatread([ dn_bsdat fn_bsdat ]);

%%% stream function/coil surface mesh
[listNode, listTriangle, tri] = importMeshWavefront('./mesh/meshGy_Example_200626.obj');


coil.listNode = listNode;
coil.listTriangle = listTriangle;
coil.tri = tri;
% 
% clear coil;
% load coil_190402_v02.mat;

%% generate target field
SZ3d = [43 43 43];


ROImsk_s1 = circshift(ROImsk_or_vol,[2 0 0]);
ROImsk_s2 = circshift(ROImsk_or_vol,[4 0 0]);

ROImsk = ROImsk_or_vol | ROImsk_s1 | ROImsk_s2;
ROImsk_er = imerode(ROImsk,ones([3 3 3]));
ROImsk_bnd = ROImsk & (~ROImsk_er);

ROImsk_tar = ROImsk_bnd;


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

xROI = x3tar(ROImsk_tar);
yROI = y3tar(ROImsk_tar);
zROI = z3tar(ROImsk_tar);

rk = cat(1, xROI(:)', yROI(:)', zROI(:)');
rk = rk';
%% Then we habe to calculate the field in a given direction


targetCoil = 'dBzdy';

Btarg = yROI * 0.005*0.1;

coil.btarget = Btarg(:);
coil.targetCoil = targetCoil;


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
reg=1.0*10^-6;
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

save('./precomp/OptSetupGy_Example_200626.mat','Btarg','ROImsk_tar','x3tar','y3tar','z3tar','optimizationType','reg','calculateR','calculateL','calculateLwp','rk','targetCoil');

save('./precomp/coil_OptSetupGy_Example_200626.mat','coil');