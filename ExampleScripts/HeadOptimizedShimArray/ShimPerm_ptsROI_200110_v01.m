%% define dimensions, points, etc;

sp_bk = 0.01;   % 1cm spacing between shim locations
Xc_max = 0.170;  % cylinder goes to 170mm along x

Rsph = 0.15;     % spherical design surface for shims has R=15cm

%% cylinder section

yc1d = 0.001* [  -16.99,  -52.63,  -80.71,  -98.05,  -98.05, -114.23, -119.62, ...
                -119.62, -114.23,  -98.05,  -98.05,  -80.71,  -52.63,  -16.99, ...
                  16.99,   52.63,   80.71,   98.05,   98.05,  114.23,  119.62, ...
                 119.62,  114.23,   98.05,   98.05,   80.71,   52.63,   16.99 ];
             
zc1d = 0.001* [  139.96,  130.79,  115.60,   96.71,   74.71,   48.44,   19.78, ...
                 -19.78,  -48.44,  -74.71,  -96.71, -115.60, -130.79, -139.96, ...
                -139.96, -130.79, -115.60,  -96.71,  -74.71,  -48.44,  -19.78, ...
                  19.78,   48.44,   74.71,   96.71,  115.60,  130.79,  139.96 ];
              
Nc_yz = numel(yc1d);
              
tc1d = atan2( yc1d, zc1d );
rc1d_yz = sqrt( yc1d.^2 + zc1d.^2 );
              
Nxc = floor( Xc_max / sp_bk );
sp_xc = Xc_max / Nxc;
xc1d = (Xc_max-sp_xc/2):(-1*sp_xc):0;

[ xc2d, yc2d ] = ndgrid( xc1d(:), yc1d(:) );
[ xc2d, zc2d ] = ndgrid( xc1d(:), zc1d(:) );
[ xc2d, tc2d ] = ndgrid( xc1d(:), tc1d(:) );

rc2d_yz = sqrt( yc2d.^2 + zc2d.^2 );
rc2d_xyz = sqrt( yc2d.^2 + zc2d.^2 + xc2d.^2 );

mask_sph = rc2d_xyz<Rsph;

rc2d_yz( mask_sph ) = sqrt( Rsph^2 - xc2d( mask_sph ).^2);
yc2d( mask_sph ) = rc2d_yz( mask_sph ) .* sin( tc2d( mask_sph ));
zc2d( mask_sph ) = rc2d_yz( mask_sph ) .* cos( tc2d( mask_sph ));

figure(11); hold off; plot3( xc2d(:), yc2d(:), zc2d(:), 'o'); grid on; grid minor; axis equal

%%% normalized Halbach magnetization vector

mxc2d = zeros(size(xc2d));
myc2d = sin(2*tc2d);
mzc2d = cos(2*tc2d);

%% sphere section
% use:
%   - tc1d
%   - Rsph
%   - phi = [-pi/2, pi/2]
ts1d = tc1d(1:Nc_yz/2);
rs1d_xy = abs(sin(ts1d))*Rsph;
Nc_p_all = floor( pi*rs1d_xy / sp_bk );

ts_all = zeros(sum(Nc_p_all(:)), 1);
ps_all = zeros(sum(Nc_p_all(:)), 1);

itmp = 1;
for irr = 1:Nc_yz/2
    Ntmp = Nc_p_all(irr);
    ts_all( itmp:(Ntmp+itmp-1) ) = ts1d(irr);
    delp = pi/Ntmp;
    ptmp_all = (-pi/2+delp/2):(delp):(pi/2);
    ps_all( itmp:(Ntmp+itmp-1) ) = ptmp_all;
    zs_all( itmp:(Ntmp+itmp-1) ) = Rsph*cos(ts1d(irr));
    itmp = itmp+Ntmp;
end

xs_all =  1*cos(ps_all).*sin(ts_all)*Rsph;
ys_all =  1*sin(ps_all).*sin(ts_all)*Rsph;

figure(11); hold on; plot3( xs_all(:), ys_all(:), zs_all(:), 'o');

%%% normalized Halbach magnetization vector

mxs2d = sin(2*ts_all).*cos(ps_all);
mys2d = sin(2*ts_all).*sin(ps_all);
mzs2d = cos(2*ts_all);


% pp_ba = cat(1, pc_all(:), ps_all(:));
% tt_ba = cat(1, tc_all(:), ts_all(:));
% ss_ba = zeros(size(pp_ba));

%% unite coordinates

xx_a = cat(1, xc2d(:), xs_all(:));
yy_a = cat(1, yc2d(:), ys_all(:));
zz_a = cat(1, zc2d(:), zs_all(:));

Nmag = numel(xx_a);

%% create initial guess halbach magnetization vectors (normalized)

mx_ba = cat(1, mxc2d(:), mxs2d(:));
my_ba = cat(1, myc2d(:), mys2d(:));
mz_ba = cat(1, mzc2d(:), mzs2d(:));

figure(12); quiver3( xx_a, yy_a, zz_a, mx_ba, my_ba, mz_ba ); axis equal

%% get target ROI information

%%% Use magnet geometry information from original magnet array design
fn_ROI = '../../ExampleScripts/HeadOptimizedMagnetArray/precomp/MagROI_190528_v01.mat';

load( fn_ROI, 'x3tar','y3tar','z3tar','ROImsk_or','ROImsk_or_vol');
ROImsk_or = ROImsk_or | flip(ROImsk_or,2);
ROImsk_or_vol = ROImsk_or_vol | flip(ROImsk_or_vol,2);

dispimgs(double(ROImsk_or),'gray','sag');
Ntar = sum(ROImsk_or(:));   
Ntar_vol = sum(ROImsk_or_vol(:));

% save('precomp/ShimROI_200110_v01.mat','xx_a','yy_a','zz_a','x3tar','y3tar','z3tar','ROImsk_or','ROImsk_or_vol','mx_ba','my_ba','mz_ba','Nmag','Ntar','Ntar_vol','tt_ba','pp_ba','ss_ba');
save('./precomp/ShimROI_200110_v01.mat','xx_a','yy_a','zz_a','x3tar','y3tar','z3tar','ROImsk_or','ROImsk_or_vol','mx_ba','my_ba','mz_ba','Nmag','Ntar','Ntar_vol');