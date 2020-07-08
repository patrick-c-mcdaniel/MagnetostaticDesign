function [ output ] = Comsol_buildmodel_func_example_200625( Nround, wRand, B_r, mu_r )
    % simplehalbach_8_24.m, halbach_3runglayer_build_full_model.m
    %
    % Clarissa Z Cooley 8/21/2018

    import com.comsol.model.*
    import com.comsol.model.util.*
    ModelUtil.clear
    ModelUtil.setServerBusyHandler(ServerBusyHandler(2));
    model = ModelUtil.create('Model');
    %

    %% load magnet M, geometry information

    %%% load magnet information (Mx,My,Mz) for all points
    %       variables:
    %           Mout   : vector of points arranged as: [ Mx1, My1, Mz1, Mx2, My2, Mz2, Mx3, ... ]
    load('./optresults/Mout_Run2Example.mat','Mout');

    %%% load magnet information (spatial positions) for all points
    %       variables:
    %           xx_a   : x-coordinates of magnet centers
    %           yy_a   : y-coordinates of magnet centers
    %           zz_a   : z-coordinates of magnet centers
    load('./precomp/MagROI_190528_v01.mat','xx_a','yy_a','zz_a');


    %%% load magnet design specification (including symmetrization + manual
    %%% adjustments
    %       variables:
    %           tht_adj : block orientation (theta; BS convention)
    %           phi_adj : block orientation (phi; BS convention)
    %           psi_adj : block orientation (psi; BS convention)
    %           xx_a    : block x-positions
    %           yy_a    : block y-positions
    %           zz_a    : block z-positions
    %           xsz_all : block size (x; BS convention); used for BS file gen.
    %           ysz_all : block size (y; BS convention); used for BS file gen.
    %           zsz_all : block size (z; BS convention); used for BS file gen.
    %           Mm      : block magnetization magnitudes (no scaling/adjustment)
    load('./optresults/MagStruct_OptExample_200625.mat','MagBlkSpec'    );
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


    Nmag = numel(xx_a);

    %% round block y-dimension
    Nrnd = Nround;
    delysz = 0.0254 / (Nrnd-1);
    ydisc = 0:delysz:0.0254;
    
    ysz_all_unrnd = ysz_all;

    ysz_all_rnd = RoundToSet( ysz_all, ydisc );
    
    ysz_randoff = wRand * (0.5 - rand(size(ysz_all_rnd)));
    ysz_all_rnd = ysz_all_rnd + ysz_randoff;
    ysz_all = ysz_all_rnd;

    xyz_a = [ xx_a(:)    yy_a(:)        zz_a(:) ];
    sz_a  = [ xsz_all(:) ysz_all(:)     zsz_all(:) ];
    
        strout = datestr(now,30);
    
    str_Nr  = [ '_N'   num2str(Nround,'%03u') ];
    str_wR  = [ '_w'   num2str(wRand*1000) ];
    str_Br  = [ '_Br'  num2str( round(B_r) ) 'p' num2str( round(100*(B_r-round(B_r))), '%02u' ) ];
    str_mur = [ '_mur' num2str( round(mu_r) ) 'p' num2str( round(100*(mu_r-round(mu_r))), '%02u' ) ];
    save([ './optresults/Mag_ydims_Rounded', str_Nr, str_wR, str_Br, str_mur, '_', strout, '.mat'], 'ysz_all_rnd','ysz_randoff');
    
    %% Use unrounded block sizes (uncomment to use rounded block sizes)
    
    ysz_all = ysz_all_unrnd;

    %% build geometry
    geom1 = model.geom.create('geom1', 3);


    sph1 = model.geom('geom1').feature.create('sph1', 'Sphere');
    sph1.set('r', '.11');
    sph1.set('pos', {'0' '0' '0'});
    sph1.set('axis', {'1' '0' '0'});

    sph3 = model.geom('geom1').feature.create('sph3', 'Sphere');
    sph3.set('r', '1');
    sph3.set('pos', {'0' '0' '0'});
    sph3.set('axis', {'1' '0' '0'});

    %place and rotate cubes
    blkcount = 0;  
    for iblk = 1:Nmag

        pos_xyz = [ xx_a(iblk) yy_a(iblk) zz_a(iblk) ];
        tht_bk = tht_adj(iblk);
        phi_bk = phi_adj(iblk);
        psi_bk = psi_adj(iblk);

        %%% note: swap y- and z- sizes because Biot-Savart swaps y and z when specifying block size (its a terrible internal inconsistency on
        %%% Biot-Savart's part - truly BS software
        xsz_bk = xsz_all(iblk);
        ysz_bk = ysz_all(iblk);
        zsz_bk = zsz_all(iblk);

        placecube_3ang_BScon(model, pos_xyz, tht_bk, phi_bk, psi_bk, xsz_bk, ysz_bk, zsz_bk, iblk);   %%function for cube placement

    end

    %%% build geometry
    model.geom('geom1').run;


    %% geometry domains (new) 

    %%% blocks
    del_tol = 1e-4;     % extra margin for bounding box

    blkcount2 = 1;
    domain_blk = zeros(1, Nmag);

    for iblk = 1:Nmag
        disp(['block number ', num2str(iblk)]);
        pos_xyz = [ xx_a(iblk) yy_a(iblk) zz_a(iblk) ];

        sz_bk   = sz_a(iblk,:);
        r_bnd   = del_tol + sqrt( sum( sz_bk.^2 ) )/2;

    %     domain_blk(blkcount2) = mphselectcoords(model,'geom1',pos_xyz,'domain','Radius',r_bnd);

        coordBox = [pos_xyz(1)-r_bnd pos_xyz(1)+r_bnd; ...
                    pos_xyz(2)-r_bnd pos_xyz(2)+r_bnd; ...
                    pos_xyz(3)-r_bnd pos_xyz(3)+r_bnd];
        domain_blk(blkcount2) = mphselectbox(model,'geom1',coordBox,'domain');

        blkcount2 = blkcount2+1;
    end



    pos1 = [0 0 0]; del_pos = 0.11;
    coordBox = [pos1(1)-del_pos pos1(1)+del_pos;pos1(2)-del_pos pos1(2)+del_pos;pos1(3)-del_pos pos1(3)+del_pos];
    domain_sphere2 = mphselectcoords(model,'geom1',pos1','domain','radius',del_pos+del_tol);

    domain_tot = [1:(Nmag+2)];
    domain_air = setdiff(domain_tot,domain_blk);

    %% geometry entity material definitions
    %
    % model.geom('geom1').run;
    model.selection.create('air');
    model.selection('air').set([domain_air]);
    model.selection.create('ndfeb');
    model.selection('ndfeb').set(domain_blk);
    %  
    model.material.create('mat1');
    model.material('mat1').selection.named('ndfeb');
    model.material('mat1').name('NdFeB');
    
    mu_r_str = num2str(mu_r);
    model.material('mat1').propertyGroup('def').set('relpermeability', {mu_r_str '0' '0' '0' mu_r_str '0' '0' '0' mu_r_str});
    % 
    model.material.create('mat2');
    model.material('mat2').selection.named('air');
    model.material('mat2').name('air');
    model.material('mat2').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    % 
    % model.cpl.create('maxop1', 'Maximum', 'geom1');
    % model.cpl.create('minop1', 'Minimum', 'geom1');
    % model.cpl('maxop1').selection.set([domain_air(2)]);
    % model.cpl('minop1').selection.set([domain_air(2)]);
    % 
    model.physics.create('mfnc', 'MagnetostaticsNoCurrents', 'geom1');


    %% mesh and study/physics setup



    %% add Br conditions into all blocks

    % note:         Mtht = atan2( Mxy, Mz )
    %               tht_z = Mtht;
    %               tht_zd = tht_z * 180/pi;
    %
    %               Mphi = atan2( My, Mx );
    %               phi_z = Mphi;
    %               phi_adj = 90+ phi_z * 180/pi

%     Br_m = 1.42;        % everything is N52 for now
    Br_m = B_r;        % everything is N52 for now
    
    Br_z  = Br_m * cos( pi/180 * tht_adj );
    Br_xy = Br_m * sin( pi/180 * tht_adj );

    phi_z = (phi_adj-90)*pi/180;
    Br_x  = cos(phi_z) .* Br_xy;
    Br_y  = sin(phi_z) .* Br_xy;


    Br_xyz = [Br_x(:) Br_y(:) Br_z(:)];

    addBr_block_181007(model, domain_blk(:), Br_xyz);


    %% mesh, study setup
    
    model.mesh.create('mesh1', 'geom1');
    model.mesh('mesh1').create('ftet1', 'FreeTet');
    model.mesh('mesh1').feature('ftet1').create('size1', 'Size');
    model.mesh('mesh1').feature('ftet1').create('size2', 'Size');
    model.mesh('mesh1').feature('ftet1').feature('size1').selection.geom('geom1', 3);
    model.mesh('mesh1').feature('ftet1').feature('size1').selection.all;
    model.mesh('mesh1').feature('ftet1').feature('size1').selection.set([2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328]);
    model.mesh('mesh1').feature('ftet1').feature('size1').set('custom', 'on');
    model.mesh('mesh1').feature('ftet1').feature('size1').set('hmaxactive', 'on');
    model.mesh('mesh1').feature('ftet1').feature('size1').set('hmax', '0.005');
    model.mesh('mesh1').feature('ftet1').feature('size1').set('hminactive', 'on');
    model.mesh('mesh1').feature('ftet1').feature('size1').set('hmin', '0.002');
    model.mesh('mesh1').feature('ftet1').feature('size2').selection.geom('geom1', 3);
    model.mesh('mesh1').feature('ftet1').feature('size2').selection.set([1]);
    model.mesh('mesh1').feature('ftet1').feature('size2').set('custom', 'on');
    model.mesh('mesh1').feature('ftet1').feature('size2').set('hmaxactive', 'on');
    model.mesh('mesh1').feature('ftet1').feature('size2').set('hmax', '0.05');
    model.mesh('mesh1').feature('ftet1').feature('size2').set('hminactive', 'on');
    model.mesh('mesh1').feature('ftet1').feature('size2').set('hmin', '0.003');
    model.mesh('mesh1').run;

    model.study.create('std1');
    model.study('std1').create('stat', 'Stationary');
    model.study('std1').feature('stat').activate('mfnc', true);

    model.sol.create('sol1');
    model.sol('sol1').study('std1');

    model.study('std1').feature('stat').set('notlistsolnum', 1);
    model.study('std1').feature('stat').set('notsolnum', '1');
    model.study('std1').feature('stat').set('listsolnum', 1);
    model.study('std1').feature('stat').set('solnum', '1');

    model.sol('sol1').create('st1', 'StudyStep');
    model.sol('sol1').feature('st1').set('study', 'std1');
    model.sol('sol1').feature('st1').set('studystep', 'stat');
    model.sol('sol1').create('v1', 'Variables');
    model.sol('sol1').feature('v1').set('control', 'stat');
    model.sol('sol1').create('s1', 'Stationary');
    model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
    model.sol('sol1').feature('s1').create('i1', 'Iterative');
    model.sol('sol1').feature('s1').feature('i1').set('linsolver', 'cg');
    model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('prefun', 'amg');
    model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'i1');
    model.sol('sol1').feature('s1').feature.remove('fcDef');
    model.sol('sol1').attach('std1');

    model.sol('sol1').runAll;

    
    %% output files
    
    model.result.export.create('data1', 'Data');
    model.result.export('data1').set('data', 'dset1');
    model.result.export('data1').setIndex('expr', 'mfnc.Bx', 0);
    model.result.export('data1').setIndex('expr', 'mfnc.By', 1);
    model.result.export('data1').setIndex('expr', 'mfnc.Bz', 2);
%     model.result.export('data1').set('alwaysask', 'on');
    model.result.export('data1').set('header', 'off');
    model.result.export('data1').set('location', 'grid');
    model.result.export('data1').set('gridx3', 'range(-0.147,0.007,0.147)');
    model.result.export('data1').set('gridy3', 'range(-0.147,0.007,0.147)');
    model.result.export('data1').set('gridz3', 'range(-0.147,0.007,0.147)');
    
    

    
    fn_dat = ['./COMSOLfiles/DataFiles/Dat_example_200625', str_Nr, str_wR, str_Br, str_mur, '_', strout, '_vol7m.txt' ];
    
    model.result.export('data1').set('filename', fn_dat);
    model.result.export('data1').run;
    

    mphsave(model, ['./COMSOLfiles/MPH_example_200625', str_Nr, str_wR, str_Br, str_mur, '_', strout, '.mph']);
    
%     save([ './matfiles/Mag_190529_v03', str_Nr, str_wR, str_Br, str_mur, '_', strout, '.mat'], 'ysz_all_rnd');
    
end