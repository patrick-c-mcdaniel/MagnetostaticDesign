function domain_n = placecube(model,pos1,rotang,blkcount)

disp(['position = [',num2str(pos1),']']);
posx = num2str(pos1(1));  
posy = num2str(pos1(2)) ;
posz = num2str(pos1(3));
rotang = num2str(rotang);

blktag = ['blk',num2str(blkcount+1)]; 


blk1 = model.geom('geom1').feature.create(blktag, 'Block');

model.geom('geom1').feature(blktag).set('size', {'1 [in]' '1 [in]' '1 [in]'});
model.geom('geom1').feature(blktag).set('pos', {posx posy posz});


model.geom('geom1').feature(blktag).set('base', 'center');


rottag = ['rot',num2str(blkcount+1)]; 
blkrot = model.geom('geom1').feature.create(rottag, 'Rotate');
model.geom('geom1').feature(rottag).set('axis', {'1' '0' '0'});
model.geom('geom1').feature(rottag).set('pos', {posx posy posz});
model.geom('geom1').feature(rottag).set('rot', rotang);
model.geom('geom1').feature(rottag).selection('input').set({blktag});

% del_pos = 2e-2;
% coordBox = [pos1(1)-del_pos pos1(1)+del_pos;pos1(2)-del_pos pos1(2)+del_pos;pos1(3)-del_pos pos1(3)+del_pos];
% domain_n = mphselectbox(model,'geom1',coordBox,'domain')


 %mphgeom(model)