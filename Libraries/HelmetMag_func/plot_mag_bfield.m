function [ outdum ] = plot_mag_bfield( Mout, optimvalues, state, Bout, ROImsk_or, xx_a, yy_a, zz_a )
%PLOT_MAG_BFIELD Summary of this function goes here
%   Detailed explanation goes here
    %% plot B=field
    Bout = D2d * Mout;
    Bm = sqrt( Bout(1:3:end).^2 + Bout(2:3:end).^2 + Bout(3:3:end).^2 );
    
    dispBvec(Bout, ROImsk_or,'mag','XZ',31);
    figure(31); colormap jet; hcb = colorbar; caxis([min(Bm(:)) max(Bm(:))]); set(hcb,'color',[1 1 1]); drawnow
    %% show magnetization
    
    Mout_x = Mout(1:3:end);
    Mout_y = Mout(2:3:end);
    Mout_z = Mout(3:3:end);

    figure(101); quiver3( xx_a, yy_a, zz_a, Mout_x, Mout_y, Mout_z ); axis equal
    
%     stopout = 0;
outdum = 0;
    save('Mout_tmp_save.mat', 'Mout');
end

