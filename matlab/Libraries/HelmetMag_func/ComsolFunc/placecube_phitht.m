function domain_n = placecube_phitht(model,pos1,tht_bk,phi_bk,blkcount)

disp(['position = [',num2str(pos1),']']);
posx = num2str(pos1(1));  
posy = num2str(pos1(2)) ;
posz = num2str(pos1(3));
tht_bk = num2str(tht_bk);
phi_bk = num2str(phi_bk);

blktag = ['blk',num2str(blkcount+1)]; 


blk1 = model.geom('geom1').feature.create(blktag, 'Block');

model.geom('geom1').feature(blktag).set('size', {'1 [in]' '1 [in]' '1 [in]'});
model.geom('geom1').feature(blktag).set('pos', {posx posy posz});


model.geom('geom1').feature(blktag).set('base', 'center');


rottag_x = ['rot_x_',num2str(blkcount+1)]; 
rottag_z = ['rot_z_',num2str(blkcount+1)]; 

blkrotx = model.geom('geom1').feature.create(rottag_x, 'Rotate');
model.geom('geom1').feature(rottag_x).set('axis', {'1' '0' '0'});
model.geom('geom1').feature(rottag_x).set('pos', {posx posy posz});
model.geom('geom1').feature(rottag_x).set('rot', tht_bk);
model.geom('geom1').feature(rottag_x).selection('input').set({blktag});

blkrotz = model.geom('geom1').feature.create(rottag_z, 'Rotate');
model.geom('geom1').feature(rottag_z).set('axis', {'0' '0' '1'});
model.geom('geom1').feature(rottag_z).set('pos', {posx posy posz});
model.geom('geom1').feature(rottag_z).set('rot', phi_bk);
model.geom('geom1').feature(rottag_z).selection('input').set({rottag_x});

% del_pos = 2e-2;
% coordBox = [pos1(1)-del_pos pos1(1)+del_pos;pos1(2)-del_pos pos1(2)+del_pos;pos1(3)-del_pos pos1(3)+del_pos];
% domain_n = mphselectbox(model,'geom1',coordBox,'domain')


 %mphgeom(model)