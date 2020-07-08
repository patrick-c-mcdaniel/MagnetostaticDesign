function [ MagBlkStruct ] = Mred2BSfile( fn_Mred, fn_MagROI, fn_MagBlkSpec, fn_BSout )

    %% Generate Biot-Savart file from GA result
    %
    %   V02 : finally symmetrized Mm...

    unix( ['touch ' fn_BSout] );

    fobj = fopen(fn_BSout,'w');

    str_head = ['Info {BiotSavart 4.1 data file}\n\n' ...
                '\n' ...
                '\n\n' ];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Put magnet blocks in place

    %% load magnet positions/values
    % 
    % load('../MagROI_180821_v01.mat');
    % load('/home3/patmcd/Work_Dir/Halbach/Magnet/Shim_180820/180823/Mout_20180823T214353.mat');

    %%% loads Brmat; contains all magnet material info
    %       0: air
    %       1: N42 (Br=1.32-ish)
    %       2: N52 (Br=1.42-ish)
    %
    % load('../ComsolModel/constructed_halbach2_brmat');

    %% load magnet info
    load( fn_Mred,'Mout');
    load( fn_MagROI );
    load( fn_MagBlkSpec,'psi_adj');
    
%     Mfull = Mred2Mxyz( Mout, inds_mid, inds_pos, inds_neg );

    %% load initial guess to get psi;

    % load('/home/patmcd/Documents/HelmetMagnet/Magnet_v2/180926_v01/MagSpec_181104_v01.mat','psi_adj'    );
    % 
    % load('Mout_20180926T193744','Mout');
    % load('MagROI_180925_v03');



    %%

    Mfull = Mred2Mxyz(Mout, inds_mid, inds_pos, inds_neg);

%     psi_adj = cat(1, psi_adj, zeros(30,1));

    [ sz3tps ] = Mxyz_to_sz3tps( Mfull, psi_adj*pi/180 );

    xsz_all = sz3tps(:,1);
    ysz_all = sz3tps(:,2);
    zsz_all = sz3tps(:,3);

    tht_adj = sz3tps(:,4) * 180/pi;
    phi_adj = sz3tps(:,5) * 180/pi + 90;
    psi_adj = sz3tps(:,6) * 180/pi - 90;

    phi_adj([1:7 92:98]) = 90*[ -1 1 1 1 -1 -1 1 1 -1 -1 1 1 -1 -1];

    Mx = Mfull(1:3:end);
    My = Mfull(2:3:end);
    Mxy = sqrt(Mx.^2 + My.^2);

    Mz = Mfull(3:3:end);


    Mm = sqrt(Mfull(1:3:end).^2 + Mfull(2:3:end).^2 + Mfull(3:3:end).^2);

    Mtht = atan2( Mxy, Mz );
    Mphi = atan2( My, Mx );

    % Nblk = numel( xx_a );

    tht_adj = Mtht * 180/pi;
    phi_z = Mphi;
    phi_adj = 90+ phi_z * 180/pi;



    Br_N42 = 1.32;
    Br_N52 = 1.42;

    Nblk = numel(xsz_all);

    % Mm_sc = Mm;
    % Mm_max = max(Mm_sc);
    % 
    % xsz_all = xsz * Mm_sc / Mm_max;
    % ysz_all = ysz * ones(size(xsz_all));
    % zsz_all = zsz * ones(size(xsz_all));


    %%
    % MUZ  = 4*pi*1e-7;               % SI units
    % 
    % Mmax_ID_1 = 2 * (0.0254)^3 * 1*1*1 / MUZ;    % SI units
    % Mmax_ID_2 = 2 * (0.0254)^3 * 0.5*0.25*0.5 / MUZ;    % SI units
    % 
    % Mout_sc = Mout / Mm_max * 1.42 * (0.0254)^3 * 1*1*1 / MUZ;
    % save('Mout_sc_181104.mat','Mout_sc');

    %%


    % % adjust psi
    % psi_adj = zeros(size(phi_zd));
    % 
    % i_ps_ctr1 = [ 1:7 ];
    % i_ps_ctr2 = [    92          93          94          95          96          97          98          ]; % numerical indices
    % 
    % a_ps_ctr = [    67         -55         -87          75          87         -85          90          ]; % angles [deg]
    % 
    % i_ps1 = [   11          12          13          14                                              ...
    %             19          20          21                                                          ...
    %             27          28                                                                      ...
    %            103         104         105                                                          ...
    %            110         111         112                                                          ...
    %             26          33                                                                      ...
    %              8         183                                                                      ];
    %        
    % i_ps2 = [  179         180         181         182                                              ...
    %            173         174         175                                                          ...
    %            167         168                                                                      ...
    %             89          90          91                                                          ...
    %             82          83          84                                                          ...
    %            159         166                                                                      ...
    %            176         185                                                                      ];
    %         
    % a_ps1 = [   20         -40         -30          30                                              ...
    %            -15         -10          18                                                          ...
    %              0           5                                                                      ...
    %             25         -10         -40                                                          ...
    %             16          -4         -20                                                          ...
    %             20          20                                                                      ...
    %             10         -15                                                                      ];
    % 
    % % psi_adj(i_ps_ctr) = a_ps_ctr;
    % psi_adj(i_ps_ctr1) = 90;
    % psi_adj(i_ps_ctr2) = 90;
    % psi_adj(i_ps1) = a_ps1;
    % psi_adj(i_ps2) = -1*a_ps1;
    % 
    % % manual edit...
    % psi_adj(27) = -5;




    %% scale block sizes so they pack without intersecting

    % i_xzsz = [  181 182 174 175 13  14  20  21  6   7   ];
    % xsz_all(i_xzsz) = xsz_all(i_xzsz) * 0.8;
    % zsz_all(i_xzsz) = zsz_all(i_xzsz) / 0.8;
    % 
    % % save('MagSpec_181104_v01.mat','tht_adj','phi_adj','psi_adj','xx_a','yy_a','zz_a','xsz_all','ysz_all','zsz_all','Mm');
    % 
    ver_all = zeros(0, 3);
    fac_all = zeros(0, 3);
    nver = 0;

    %% iterate through magnets; write to BSfile
    for iblk = 1:Nblk

        xSize = xsz_all(iblk);
        ySize = ysz_all(iblk);
        zSize = zsz_all(iblk);

    %     if abs(zSize) < t_tol
    %         continue;
    %     end

        xPos = xx_a(iblk);
        yPos = yy_a(iblk);
        zPos = zz_a(iblk);

        tht  = tht_adj(iblk);
        phi  = phi_adj(iblk);
        psi  = psi_adj(iblk);


        [ fac_blk ver_blk ] = gen_facver([ xPos, yPos, zPos],[ xSize ySize zSize ], phi*pi/180, tht*pi/180, psi*pi/180);
        fac_all = cat(1, fac_all, nver+fac_blk);
        ver_all = cat(1, ver_all, ver_blk);
        nver = nver + 8;

    %     str_mat = 'matdef';

        str_mat = 'N52';

      %  disp(['iblk = ' num2str(iblk) '  phideg = ' num2str(phideg) '  thtdeg = ' num2str(thtdeg) '  psideg = ' num2str(psideg) ]);
        tmp_str = ['Cuboid {\n' ...
                   'name {Mag' num2str(iblk,'%04u') '}\n' ...
                   'color 19660 45874 45874\n' ...
                   'positionX ' num2str(xPos) '\n' ...
                   'positionY ' num2str(yPos) '\n' ...
                   'positionZ ' num2str(zPos) '\n' ...
                   'eulerPhi '   num2str(phi) '\n' ...
                   'eulerTheta ' num2str(tht) '\n' ...
                   'eulerPsi '   num2str(psi) '\n' ...
                   'currentSupply none\n' ...
                   'wireDiameter 0\n' ...
                   'winding 1\n' ...
                   'material {' str_mat '}\n' ...
                   'resolution 0.5mm\n' ...
                   'latticeDrawType {0}\n' ...
                   'maximumMultipoleOrder 5\n' ...
                   'maximumMultipoleOrderForLeaf 1\n' ...
                   'drawEdge\n' ...
                   'drawSurfaceElements\n' ...
                   'sizeX ' num2str(xSize) '\n' ...
                   'sizeY ' num2str(ySize) '\n' ...
                   'sizeZ ' num2str(zSize) '\n' ...
                   '}\n' ];
        fprintf(fobj, tmp_str);
    end

    %% magnet materials
    tmp_str = [ 'HardMaterial {\n' ...
                'name {vac}\n' ...
                'remanence 0.000000e+000\n' ...
                'remanenceDirector 0 0 1\n' ...
                'normalizeRemanenceDirector\n' ...
                '}\n\n' ];
    fprintf(fobj, tmp_str);


    tmp_str = [ 'HardMaterial {\n' ...
                'name {N42}\n' ...
                'remanence ' num2str(Br_N42) '\n' ...
                'remanenceDirector 0 0 1\n' ...
                'normalizeRemanenceDirector\n' ...
                '}\n\n' ];
    fprintf(fobj, tmp_str);

    tmp_str = [ 'HardMaterial {\n' ...
                'name {N52}\n' ...
                'remanence ' num2str(Br_N52) '\n' ...
                'remanenceDirector 0 0 1\n' ...
                'normalizeRemanenceDirector\n' ...
                '}\n\n' ];
    fprintf(fobj, tmp_str);
    %     if strcmp(str_mat, 'MatTestPos')
    %         matimg(ilp) = 1;
    %     end
    fclose(fobj);

    figure(11); patch('Faces',fac_all,'Vertices',ver_all,'FaceColor',[1 0.7 0.7])
    axis equal
    %% put block specification variables into MagBlkStruct
    
    MagBlkStruct.tht_adj    = tht_adj;  % note: angles in degrees
    MagBlkStruct.phi_adj    = phi_adj;  % note: angles in degrees
    MagBlkStruct.psi_adj    = psi_adj;  % note: angles in degrees
    MagBlkStruct.xx_a       = xx_a;
    MagBlkStruct.yy_a       = yy_a;
    MagBlkStruct.zz_a       = zz_a;
    MagBlkStruct.xsz_all    = xsz_all;
    MagBlkStruct.ysz_all    = ysz_all;
    MagBlkStruct.zsz_all    = zsz_all;
    MagBlkStruct.Mm         = Mm;
    
%     stlwrite('Mag_v2_190521_v01.stl', fac_all, ver_all);

end