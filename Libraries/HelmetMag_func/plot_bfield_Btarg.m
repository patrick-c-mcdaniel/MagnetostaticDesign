function [ outdum ] = plot_bfield_Btarg( Mout, optimvalues, state, Bt, ROImsk_or, cm_rb, figoff )
%PLOT_MAG_BFIELD Summary of this function goes here
%   Detailed explanation goes here
%     if ishandle(100+figoff)
%         close((100+figoff));
%     end
%     if ishandle(30+figoff)
%         close((30+figoff));
%     end
    
    if nargin<7
        figoff=1;
    end
    %% plot B=field
%     Bout = D2d * Mout;

    Bmat = zeros(size(ROImsk_or));
    Bmat(ROImsk_or) = Bt;
    
    mosaic(permute( Bmat, [1 3 2]) , 6, 8, 30+figoff); figure(30+figoff); colormap(cm_rb); colorbar; caxis(mean(Bt(:))+[-5e-4 5e-4]); drawnow
%     %% show magnetization
%     
%     Mout_x = Mout(1:3:end);
%     Mout_y = Mout(2:3:end);
%     Mout_z = Mout(3:3:end);
% 
%     figure(100+figoff); hold off; quiver3( xx_a, yy_a, zz_a, Mout_x, Mout_y, Mout_z ); axis equal
%     drawnow;
% %     stopout = 0;
% outdum = 0;
    save('Mout_tmp_save.mat', 'Mout');
end

