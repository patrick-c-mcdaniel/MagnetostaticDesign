function domain_n = placecube_3ang_BScon(model,pos1,tht_bk,phi_bk,psi_bk, xsz_bk, ysz_bk, zsz_bk, blkcount)

disp(['position = [',num2str(pos1),']']);
posx = num2str(pos1(1));  
posy = num2str(pos1(2)) ;
posz = num2str(pos1(3));
tht_bk = num2str(tht_bk);
phi_bk = num2str(phi_bk);
psi_bk = num2str(psi_bk);

xsz_str = num2str(xsz_bk);
ysz_str = num2str(ysz_bk);
zsz_str = num2str(zsz_bk);

blktag = ['blk',num2str(blkcount+1)]; 


blk1 = model.geom('geom1').feature.create(blktag, 'Block');

model.geom('geom1').feature(blktag).set('size', {[xsz_str ' [m]'] [ysz_str ' [m]'] [zsz_str ' [m]']});
model.geom('geom1').feature(blktag).set('pos', {posx posy posz});


model.geom('geom1').feature(blktag).set('base', 'center');

rottag_zpp = ['rot_zpp_',num2str(blkcount+1)];  % psi rotation 
rottag_xp  = ['rot_xp_',num2str(blkcount+1)];   % theta rotation
rottag_z   = ['rot_z_',num2str(blkcount+1)];    % phi rotation

blkrot_zpp = model.geom('geom1').feature.create(rottag_zpp, 'Rotate');
model.geom('geom1').feature(rottag_zpp).set('axis', {'0' '0' '1'});
model.geom('geom1').feature(rottag_zpp).set('pos', {posx posy posz});
model.geom('geom1').feature(rottag_zpp).set('rot', psi_bk);
model.geom('geom1').feature(rottag_zpp).selection('input').set({blktag});

blkrot_xp = model.geom('geom1').feature.create(rottag_xp, 'Rotate');
model.geom('geom1').feature(rottag_xp).set('axis', {'1' '0' '0'});
model.geom('geom1').feature(rottag_xp).set('pos', {posx posy posz});
model.geom('geom1').feature(rottag_xp).set('rot', tht_bk);
model.geom('geom1').feature(rottag_xp).selection('input').set({rottag_zpp});

blkrot_z = model.geom('geom1').feature.create(rottag_z, 'Rotate');
model.geom('geom1').feature(rottag_z).set('axis', {'0' '0' '1'});
model.geom('geom1').feature(rottag_z).set('pos', {posx posy posz});
model.geom('geom1').feature(rottag_z).set('rot', phi_bk);
model.geom('geom1').feature(rottag_z).selection('input').set({rottag_xp});

% del_pos = 2e-2;
% coordBox = [pos1(1)-del_pos pos1(1)+del_pos;pos1(2)-del_pos pos1(2)+del_pos;pos1(3)-del_pos pos1(3)+del_pos];
% domain_n = mphselectbox(model,'geom1',coordBox,'domain')


 %mphgeom(model)