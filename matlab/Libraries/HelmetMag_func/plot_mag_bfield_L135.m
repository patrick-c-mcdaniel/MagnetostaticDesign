function [ outdum ] = plot_mag_bfield( Mout, optimvalues, state, Bout, ROImsk_or, xx_a, yy_a, zz_a,cm_rb, figoff )
%PLOT_MAG_BFIELD Summary of this function goes here
%   Detailed explanation goes here
%     if ishandle(100+figoff)
%         close((100+figoff));
%     end
%     if ishandle(30+figoff)
%         close((30+figoff));
%     end
    
    if nargin<10
        figoff=1;
    end
    %% plot B=field
%     Bout = D2d * Mout;
    Bm = sqrt( Bout(1:3:end).^2 + Bout(2:3:end).^2 + Bout(3:3:end).^2 );
    
    dispBvec(Bout, ROImsk_or,'mag','XZ',31);
    figure(30+figoff); hold off; colormap(cm_rb); hcb = colorbar; caxis(mean(Bm(:)) + [-5e-4 5e-4]); set(hcb,'color',[1 1 1]); drawnow
    %% show magnetization
    
    Mout_x = Mout(1:3:end);
    Mout_y = Mout(2:3:end);
    Mout_z = Mout(3:3:end);

    figure(100+figoff); hold off; quiver3( xx_a, yy_a, zz_a, Mout_x, Mout_y, Mout_z ); axis equal
    drawnow;
%     stopout = 0;
outdum = 0;
    save('Mout_tmp_save.mat', 'Mout');
end

