function [ D2d ] = proc_Dmatrix_vol_L1( fn_MagROI, fn_D1save )
    load(fn_MagROI);

    xROI = x3tar( ROImsk_or_vol );
    yROI = y3tar( ROImsk_or_vol );
    zROI = z3tar( ROImsk_or_vol );

    Nmag = numel(xx_a);
    Ntar = numel(xROI);

    MUZ = 4*pi*1e-7;   % vacuum permeability

    %% Create forward calculation dipole matrix
    %
    % [xROI yROI zROI] and [xx_a' yy_a' zz_a'] contain ROI and magnet
    %    coordinates, respectively. 
    % Create 

    D4d = zeros(3, 3, Ntar, Nmag);
    D2d = zeros(3*Ntar, 3*Nmag);

    % compute "outer difference" between sets of points:

    x_mt = repmat( xROI, [1 Nmag] ) - repmat( xx_a',[Ntar 1] );
    y_mt = repmat( yROI, [1 Nmag] ) - repmat( yy_a',[Ntar 1] );
    z_mt = repmat( zROI, [1 Nmag] ) - repmat( zz_a',[Ntar 1] );
    r_mt = sqrt( x_mt.^2 + y_mt.^2 + z_mt.^2 );

    % fill 4D Dipole matrix
    D4d(1, 1, :, :) = MUZ/(4*pi) * reshape( -1./(r_mt.^3) + 3*(x_mt.*x_mt)./(r_mt.^5), [1 1 Ntar Nmag] );
    D4d(1, 2, :, :) = MUZ/(4*pi) * reshape(                 3*(x_mt.*y_mt)./(r_mt.^5), [1 1 Ntar Nmag] );
    D4d(1, 3, :, :) = MUZ/(4*pi) * reshape(                 3*(x_mt.*z_mt)./(r_mt.^5), [1 1 Ntar Nmag] );

    D4d(2, 1, :, :) = MUZ/(4*pi) * reshape(                 3*(y_mt.*x_mt)./(r_mt.^5), [1 1 Ntar Nmag] );
    D4d(2, 2, :, :) = MUZ/(4*pi) * reshape( -1./(r_mt.^3) + 3*(y_mt.*y_mt)./(r_mt.^5), [1 1 Ntar Nmag] );
    D4d(2, 3, :, :) = MUZ/(4*pi) * reshape(                 3*(y_mt.*z_mt)./(r_mt.^5), [1 1 Ntar Nmag] );

    D4d(3, 1, :, :) = MUZ/(4*pi) * reshape(                 3*(z_mt.*x_mt)./(r_mt.^5), [1 1 Ntar Nmag] );
    D4d(3, 2, :, :) = MUZ/(4*pi) * reshape(                 3*(z_mt.*y_mt)./(r_mt.^5), [1 1 Ntar Nmag] );
    D4d(3, 3, :, :) = MUZ/(4*pi) * reshape( -1./(r_mt.^3) + 3*(z_mt.*z_mt)./(r_mt.^5), [1 1 Ntar Nmag] );

    % turn into 2D Dipole matrix (see slides to understand matrix structure)

    D2d(1:3:end, 1:3:end) = D4d(1, 1, :, :);
    D2d(1:3:end, 2:3:end) = D4d(1, 2, :, :);
    D2d(1:3:end, 3:3:end) = D4d(1, 3, :, :); 

    D2d(2:3:end, 1:3:end) = D4d(2, 1, :, :);
    D2d(2:3:end, 2:3:end) = D4d(2, 2, :, :);
    D2d(2:3:end, 3:3:end) = D4d(2, 3, :, :); 

    D2d(3:3:end, 1:3:end) = D4d(3, 1, :, :);
    D2d(3:3:end, 2:3:end) = D4d(3, 2, :, :);
    D2d(3:3:end, 3:3:end) = D4d(3, 3, :, :); 

    if isempty(fn_D1save)
        return
    end
    save(fn_D1save,'D2d','D4d','Nmag','Ntar');
end