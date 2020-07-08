function [ ver_all, fac_all ] = genannumeshz( zcoord, Rmin, Rmax, Nr, Nc, dispopt )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    if nargin<6
        dispopt=0;
    end
    r1d = linspace(Rmin, Rmax, Nr+1);
    p1d = linspace(0, 2*pi, Nc+1);
    p1d = p1d(1:end-1);
    delp = p1d(2)-p1d(1);
    
    z1d = zcoord;
    
    x_all = zeros(1,0);
    y_all = zeros(1,0);
    z_all = zeros(1,0);
    
    for irr = 1:(Nr+1)
        rad = r1d(irr);
        x1d = rad*cos(p1d+(irr-1)*delp/2);
        y1d = rad*sin(p1d+(irr-1)*delp/2);

        x_ll = repmat( x1d(:), [1, 1] );       x_all = cat(1, x_all(:), x_ll(:));
        y_ll = repmat( y1d(:), [1, 1] );       y_all = cat(1, y_all(:), y_ll(:));
        z_ll = repmat( z1d, [1, Nc  ] )'; z_all = cat(1, z_all(:), z_ll(:));
    end
    
    ver_all = zeros(numel(x_all),3);
    ver_all(:,1) = x_all;
    ver_all(:,2) = y_all;
    ver_all(:,3) = z_all;

    %% create facets
    
    fac_all = zeros(0, 3);
    
    for irr = 1:Nr
        
        fac_tmp1 = (irr-1)*Nc+[ (1:Nc)'               , mod((1:Nc)',Nc)+1  , Nc+(1:Nc)'        ];
        fac_tmp2 = (irr-1)*Nc+[ mod((1:Nc)',Nc)+1+Nc  , Nc+(1:Nc)'         , mod((1:Nc)',Nc)+1 ];
        
        fac_tmp = cat(1, fac_tmp1, fac_tmp2);
        fac_all = cat(1, fac_all, fac_tmp);
        
    end
    
    if Rmin > Rmax
        fac_all_adj = fac_all;
        fac_all_adj(:,2) = fac_all(:,3);
        fac_all_adj(:,3) = fac_all(:,2);
        fac_all = fac_all_adj;
    end
    
    if dispopt
        patchsexy( 'Faces',fac_all,'Vertices',ver_all);
    end
end

