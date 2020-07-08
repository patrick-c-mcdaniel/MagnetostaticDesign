function [ ver_all, fac_all, Rot_eul_all ] = Mspec2blocks_plug_booster( MagBlkSpec, tol, hplug, facnums, fignum )    
    %%%% generate block model of magnet block set; store as mesh surface (each block
    %    is N=12 triangles); plot mesh;

    %% block model of initial guess magnetization

%     load( fn_Mout,'Mout');
%     load( fn_MagROI );
% 
%     Mfull = Mred2Mxyz( Mout, inds_mid, inds_pos, inds_neg );
% 
% %     load('/home/patmcd/Documents/HelmetMagnet/Magnet_v2/180926_v01/MagSpec_181104_v01.mat','tht_adj','phi_adj','psi_adj','Mm','xsz_all','ysz_all','zsz_all'    );
%     load( fn_MagBlkSpec, 'psi_adj'  );
% 
% 
%     psi_adj = cat(1, psi_adj, zeros(30,1));
% 
%     [ sz3tps ] = Mxyz_to_sz3tps( Mfull, psi_adj*pi/180 );
% 
%     xsz_all = sz3tps(:,1);
%     ysz_all = sz3tps(:,2);
%     zsz_all = sz3tps(:,3);
% 
%     tht_adj = sz3tps(:,4) * 180/pi;
%     phi_adj = sz3tps(:,5) * 180/pi + 90;
%     psi_adj = sz3tps(:,6) * 180/pi - 90;
% 
%     phi_adj([1:7 92:98]) = 90*[ -1 1 1 1 -1 -1 1 1 -1 -1 1 1 -1 -1];
if nargin<5
    fignum=11;
end

    tht_adj     = MagBlkSpec.tht_adj;
    phi_adj     = MagBlkSpec.phi_adj;
    psi_adj     = MagBlkSpec.psi_adj;
    xx_a        = MagBlkSpec.xx_a;
    yy_a        = MagBlkSpec.yy_a;
    zz_a        = MagBlkSpec.zz_a;
    xsz_all     = MagBlkSpec.xsz_all;
    ysz_all     = MagBlkSpec.ysz_all;
    zsz_all     = MagBlkSpec.zsz_all;
    Mm          = MagBlkSpec.Mm;

    %% generate block model; display

    ver_all = zeros(0, 3);
    fac_all = zeros(0, 3);
    nver = 0;
Rot_eul_all = zeros(numel(xx_a), 3, 3);

    for iblk = 1:numel(xx_a)    

        xSize = xsz_all(iblk) + 2*tol;
        ySize = ysz_all(iblk) + 2*tol;
        zSize = zsz_all(iblk) + 2*tol;

    %     if abs(zSize) < t_tol
    %         continue;
    %     end

        xPos = xx_a(iblk);
        yPos = yy_a(iblk);
        zPos = zz_a(iblk);

        tht  = tht_adj(iblk);
        phi  = phi_adj(iblk);
        psi  = psi_adj(iblk);


        [ fac_blk ver_blk, Rot_eul ] = gen_facver_plug([ xPos, yPos, zPos],[ xSize ySize zSize ], phi*pi/180, tht*pi/180, psi*pi/180, facnums(iblk), hplug);
        fac_all = cat(1, fac_all, nver+fac_blk);
        ver_all = cat(1, ver_all, ver_blk);
        Rot_eul_all(iblk, :, :) = Rot_eul;
        nver = nver + 8;
    end
    figure(fignum); patch('Faces',fac_all,'Vertices',ver_all,'FaceColor',[0.6 0.6 0.6]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis image

end