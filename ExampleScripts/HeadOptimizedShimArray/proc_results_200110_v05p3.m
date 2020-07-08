Nrun = 3;
load('./precomp/D1comp_200110_v01.mat');
load('./precomp/ShimROI_200110_v01.mat');

MUz = 4*pi*1e-7;
Br = 1.32;
BrN45 = 1.38;

Mdel = 1*1/MUz*Br*(1/8)^3*(0.0254^3);
M_N45_2x = 1/MUz*BrN45*(1/4*1/8*1/8)*(0.0254^3);

Mopt = [ 0:2 ] * M_N45_2x;

B0_mm_all = zeros(Nrun, 1);

Npop = 200;
initpop = zeros(Npop, 890);

Mred_all = zeros(Nrun, 890);

for irr = 1:Nrun
    
    tmp = ls([ './optresults/Mout_' num2str(irr, '%03u') '*' ]);
    tmp2 = load(tmp(1:end-1));
    
    Mred = tmp2.Mout;
    
    B_bnd = D1comp*Mopt(Mred+1)' + Btarg;
    
    B0_mm_all(irr) = max(B_bnd(:)) - min(B_bnd(:));
    
    Mred_all( irr,: ) = Mred;
    
end

[ B0_mm_sort, ind_B0mm ] = sort( B0_mm_all,'ascend');

for irr = ind_B0mm(1)
    
    tmp = ls([ './optresults/Mout_' num2str(irr, '%03u') '*' ]);
    tmp2 = load(tmp(1:end-1));
    
    Mred = tmp2.Mout;
    
    B_vol = D1comp_vol*Mopt(Mred+1)' + Btarg_vol;
%     
%     B0_mm_all(irr) = max(B_bnd(:)) - min(B_bnd(:));
%     
%     Mred_all( irr,: ) = Mred;
    
end

% initpop(1:50,:) = Mred_all(ind_B0mm(1:50),:);
% 
% initpop(51:100,:) = initpop(1:50,:); % + round(-1+2*rand(size(initpop(1:50,:))));
% initpop(101:150,:) = initpop(1:50,:); % + round(-1+2*rand(size(initpop(1:50,:))));
% initpop(151:200,:) = initpop(1:50,:); % + round(-1+2*rand(size(initpop(1:50,:))));
% 
% % initpop(initpop>4) = 4;
% initpop(initpop<0) = 0;
% 
% save('initpop_190626_v02.mat','initpop');

%% plot best shim result

Mxyz_opt = Mred_all( ind_B0mm(1), :)';

Mx_opt = mx_ba.*Mxyz_opt;
My_opt = my_ba.*Mxyz_opt;
Mz_opt = mz_ba.*Mxyz_opt;

figure(120); quiver3( xx_a, yy_a, zz_a, Mx_opt, My_opt, Mz_opt ); axis equal;

figure(121); histogram( Mxyz_opt, -5.5 + [1:1:10]); grid on; grid minor;

muBshim = mean(B_vol(:));

dispBmvec( B_vol,     ROImsk_or_vol, 'mag', 'XZ', 122 ); colormap(cm_rb); caxis(muBshim+[-5e-4 5e-4]);
% dispBmvec( B_vol,     ones(size(ROImsk_or_vol)), 'mag', 'XZ', 123 ); colormap(cm_rb); caxis(muBshim+[-5e-4 5e-4]);

Bm_ROI_ext = B_vol;
Bm_ROI_ext_mat = zeros(43,43,43);
Bm_ROI_ext_mat(ROImsk_or_vol) = B_vol;
Bx_ROI = zeros(size(B_vol));
By_ROI = zeros(size(B_vol));
Bz_ROI = B_vol;

save('./optresults/Bxyz_ShimExampleBest_200626.mat','Bm_ROI_ext','Bm_ROI_ext_mat','Bx_ROI','By_ROI','Bz_ROI');

%% create BS shim blocks
Mred_simp = Mred;
Mred_simp(Mred_simp==1) = 0;

    B_vol_simp = D1comp_vol*Mopt(Mred_simp+1)' + Btarg_vol;
    B0_mm_simp = max(B_vol_simp(:)) - min(B_vol_simp(:));



Mall = Mopt(Mred_simp+1)';
Nblk = numel(Mall);

Mfull = mag2vec(Mx_opt, My_opt, Mz_opt);

xsz = 1/4*0.0254*ones(Nblk,1);
ysz = 1/4*0.0254*ones(Nblk,1);
zsz = Mall / Mopt(end) * 1/8*0.0254;

pps = atan2( my_ba, mx_ba );
tts = atan2( sqrt(mx_ba.^2 + my_ba.^2), mz_ba );

sss = zeros(size(pps));


sz3tps = [ xsz(:) ysz(:) zsz(:) tts(:) pps(:) sss(:) ];
xyz_pos = [ xx_a(:), yy_a(:), zz_a(:) ];
fn_BSout = './BSfiles/BSmodel_Shim_Raw_200626.biot';

MagBlkStruct = Mfull2BSfile( Mfull(:), sz3tps, xyz_pos, fn_BSout );

%% create reduced object with no "zero" blocks
ii_nz = find( MagBlkStruct.zsz_all > 0 );

% MagBlkStruct_nz = MagBlkStruct;
MagBlkStruct_nz_xoff.Mm = MagBlkStruct.Mm(ii_nz);
MagBlkStruct_nz_xoff.xx_a = MagBlkStruct.xx_a(ii_nz);
MagBlkStruct_nz_xoff.yy_a = MagBlkStruct.yy_a(ii_nz);
MagBlkStruct_nz_xoff.zz_a = MagBlkStruct.zz_a(ii_nz);
MagBlkStruct_nz_xoff.tht_adj = MagBlkStruct.tht_adj(ii_nz);
MagBlkStruct_nz_xoff.phi_adj = MagBlkStruct.phi_adj(ii_nz);
MagBlkStruct_nz_xoff.psi_adj = MagBlkStruct.psi_adj(ii_nz);
MagBlkStruct_nz_xoff.xsz_all = MagBlkStruct.xsz_all(ii_nz);
MagBlkStruct_nz_xoff.ysz_all = MagBlkStruct.ysz_all(ii_nz);
MagBlkStruct_nz_xoff.zsz_all = MagBlkStruct.zsz_all(ii_nz);

%%% move x=165mm blocks back 1 mm
msk_xmax = MagBlkStruct_nz_xoff.xx_a > 0.164;
MagBlkStruct_nz_xoff.xx_a(msk_xmax) = 0.164;

%% create STL (with tolerance)
TOL = 0.1e-3;
hplug = 20e-3;
facnums = 4*ones(size(MagBlkStruct_nz_xoff.ysz_all));
ZCUT = 0.0935;

facnums( MagBlkStruct_nz_xoff.zz_a > ZCUT )  = 3;
facnums( MagBlkStruct_nz_xoff.zz_a <= ZCUT & MagBlkStruct_nz_xoff.zz_a >0) = 2;
facnums( MagBlkStruct_nz_xoff.zz_a < -ZCUT ) = 1;
facnums( msk_xmax ) = 5;
facnums( msk_xmax & (MagBlkStruct_nz_xoff.yy_a .* MagBlkStruct_nz_xoff.zz_a)<0 ) = 6;
fignum = 201;

 [ ver_all, fac_all ] = Mspec2blocks_plug( MagBlkStruct_nz_xoff, TOL, hplug, facnums, fignum )   ;
 
 stlwrite( './STLfiles/STLmodel_shim_200110_v05p3.stl', fac_all, ver_all );
 
 %% create BS file (for blocks; no tolerance)
 
 fn_BSout_nz = 'BSfiles/BSmodel_Shim_NoZero_200626.biot';
 
%  Mfull_nz = mag2vec( MagBlkStruct_nz.
 sz3tps_nz = [  MagBlkStruct_nz_xoff.xsz_all, ...
                MagBlkStruct_nz_xoff.ysz_all, ...
                MagBlkStruct_nz_xoff.zsz_all, ...
                MagBlkStruct_nz_xoff.tht_adj * pi/180, ...
               (MagBlkStruct_nz_xoff.phi_adj -90 ) * pi/180, ...
               (MagBlkStruct_nz_xoff.psi_adj +90 ) * pi/180 ];

 xyz_pos_nz = [ MagBlkStruct_nz_xoff.xx_a, ...
                MagBlkStruct_nz_xoff.yy_a, ...
                MagBlkStruct_nz_xoff.zz_a ];
            
MagBlkStruct_tmp = Mfull2BSfile( zeros(size(MagBlkStruct_nz_xoff.xx_a)), sz3tps_nz, xyz_pos_nz, fn_BSout_nz );











