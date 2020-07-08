function [  ] = plot_image_points( x1d, y1d, B1d, Bminmax, cmap, fignum, mshape )

%%% plot Bm_map_corr at arbitrary points (in color) at an XZ plane with
%%% some offset from the origin

%%% Inputs:
%
%   xx_gs, yy_gs, zz_gs :       vectors with coordinates of mapped points
%   Bm_map_corr :               vectors with measurement from each point
%   Bminmax :                   2x1 vector with min/max for color scale
%   yoff :                      y-location at which to plot field
%   ytol :                      tolerance (goes each way, so really half
%                                   the slice thickness)
%   fignum :                      figure number
%
%
%
%
%
%

    if nargin < 7
        mshape = 's';
    end
    
    Ncmap = size(cmap,1);
    Bm_min = Bminmax(1);
    Bm_max = Bminmax(2);
    
    fm_ind = round(Ncmap * (B1d - Bm_min) / (Bm_max - Bm_min));
    fm_ind(fm_ind<1) = 1;
    fm_ind(fm_ind>Ncmap) = Ncmap;
            
    figure(fignum)

    hold off;
    hbmag = scatter(x1d, y1d,64,[0 0 1],mshape,'filled'); % axis square;
%             tmp = colormap('jet');
    
    tmp = colormap(cmap);
    cdat_all = tmp(fm_ind,:);
    hbmag.CData = cdat_all;
    set(hbmag, 'MarkerEdgeColor','black');
    axis equal
    grid on; grid minor; colorbar; caxis([Bm_min Bm_max]);
    drawnow


end