function [  ] = plot_bfield_trajptsXY( xx_gs, yy_gs, zz_gs, Bm_map_corr, Bminmax , zoff, ztol, fignum )

%%% plot Bm_map_corr at arbitrary points (in color) at an XY plane with
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

    %% plot field map w/ traj matrix coordinates
    %%% y = 0
    xx_vis_all = zeros(1,0);
    yy_vis_all = zeros(1,0);
    ii_vis_all = zeros(1,0);
    iy0=1;

    for ipt = 1:numel(zz_gs)
        if abs(zz_gs(ipt) - zoff)<=ztol

            xx_vis_all(iy0) = xx_gs(ipt);
            yy_vis_all(iy0) = yy_gs(ipt);
            ii_vis_all(iy0) = ipt;
            iy0 = iy0 + 1;


            field_mag_T = Bm_map_corr(ii_vis_all);

            Bm_max = Bminmax(2);
            Bm_min = Bminmax(1);
            fm_ind = round(64 * (field_mag_T - Bm_min) / (Bm_max - Bm_min));
            fm_ind(fm_ind<1) = 1;
            fm_ind(fm_ind>64) = 64;


        end
    end

            figure(fignum);
    %         clf(1);
        %     subplot(2,1,1); set(gca,'xdir','reverse');

            hold off;
            hbmag = scatter(xx_vis_all, yy_vis_all,'filled'); % axis square;
            tmp = colormap('jet');
            cdat_all = tmp(fm_ind,:);
            hbmag.CData = cdat_all;
            axis equal
            grid on; grid minor; colorbar; caxis([Bm_min Bm_max]);
            drawnow




end