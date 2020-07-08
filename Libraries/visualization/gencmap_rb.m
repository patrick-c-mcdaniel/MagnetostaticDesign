
cmjet = colormap('jet');

cm2 = cat(1, cmjet(1:24,:), cmjet(41:end,:));

cm2(9:24,1) = cm2(9:24,2);
cm2(24:(end-9),3) = cm2(24:(end-9),2);

cm_rb1 = cm2;
cm_rb1 = cat(1, [0 0 0.53], cm_rb1);

rvals = interp1(cm_rb1(:,1), [1:0.5:size(cm_rb1,1)]);
gvals = interp1(cm_rb1(:,2), [1:0.5:size(cm_rb1,1)]);
bvals = interp1(cm_rb1(:,3), [1:0.5:size(cm_rb1,1)]);

% clear cm_rb1;
cm_rb = cat(2, rvals(:), gvals(:), bvals(:));


save('cm_rb.mat','cm_rb');