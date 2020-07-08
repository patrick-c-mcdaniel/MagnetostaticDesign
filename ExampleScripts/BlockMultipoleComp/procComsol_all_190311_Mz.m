Nfn = 20;
CHI = 0.05;
MUZ = 4*pi*1e-7;

lall = 0:5;
Nl=numel(lall);

coeffs = zeros(Nfn+1, Nl, 2*Nl-1);

for ifn =  1:Nfn
    disp([ 'ifn = ' num2str(ifn) ' / ' num2str(Nfn)]);
    %% load file
    
    str_fn = [ './Data/Block_y' num2str(ifn, '%02u') '_Bxyz_190311_Mz.txt' ];
    tmp = dlmread(str_fn,'',9,0);
    Bz = tmp(:,4);
    
    zdim = numel(Bz)/200/206;
    
    SZ3d   = [ 200 zdim 206 ];
    SZ3d_c = [ 200 zdim 204 ];
    
    Bz_rs = reshape(Bz, SZ3d);
    Bz_rs = 1/8*( Bz_rs + ...
                  flip( Bz_rs, 1)           + flip( Bz_rs, 2)           + flip( Bz_rs, 3)           + ...
                  flip(flip( Bz_rs, 1),2)   + flip(flip( Bz_rs, 1),3)   + flip(flip( Bz_rs, 2),3)   + ...
                  flip(flip(flip( Bz_rs, 1), 2), 3) );
    
    Bz_rs = smooth3(Bz_rs);
%     figure(11); imagesc(squeeze(Bz_rs(76,:,:))); axis image; colormap jet
%     figure(12); imagesc(squeeze(Bz_rs(:,76,:))); axis image; colormap jet
%     figure(13); imagesc(squeeze(Bz_rs(:,:,(zdim+1)/2))); axis image; colormap jet
%     
    



    %% generate 


%     xr3 = xx3;
%     yr3 = cos( deg2rad( ANG ))*yy3 - sin( deg2rad( ANG ))*zz3;
%     zr3 = cos( deg2rad( ANG ))*zz3 + sin( deg2rad( ANG ))*yy3;
    
    x3d = tmp(:,1);
    y3d = tmp(:,2);
    z3d = tmp(:,3);
    x3d_rs = reshape(x3d, SZ3d); % x3d_rs = 
    y3d_rs = reshape(y3d, SZ3d);
    z3d_rs = reshape(z3d, SZ3d);
    
    x3d_c = x3d_rs(:,:,2:(end-1));
    y3d_c = y3d_rs(:,:,2:(end-1));
    z3d_c = z3d_rs(:,:,2:(end-1));
    clear x3d y3d z3d;
    x3d = x3d_c(:);
    y3d = y3d_c(:);
    z3d = z3d_c(:);
    

    x1d = squeeze(x3d_c(:,1,1));
    y1d = squeeze(y3d_c(1,:,1));
    z1d = squeeze(z3d_c(1,1,:));
    
    DX = x1d(2)-x1d(1);
    DY = y1d(2)-y1d(1);
    DZ = z1d(2)-z1d(1);

    %%% look at z-difference;

    delz = DZ;
    dzBz_all = zeros(size(Bz_rs));
    dzBz_all(:,:,2:(end-1)) = (Bz_rs(:,:,3:end) - Bz_rs(:,:,1:(end-2))) / (2*delz);
    dzBz = -1*dzBz_all(:,:,2:(end-1));
    load('cm_rb.mat');

    %%% use for 3D spherical harmonic computation for integrals
    % rr3 = sqrt(xx3.^2 + yy3.^2 + zz3.^2);
    % tt3 = acos( zz3./(rr3) );
    % pp3 = atan2( yy3, xx3 );

    r3d = sqrt(x3d.^2 + y3d.^2 + z3d.^2);
    t3d = acos( z3d./(r3d) );
    p3d = atan2( y3d, x3d );




    for inl = 1:Nl

        ill = lall(inl);
        mall = (-1*ill):(ill);
        Nm = numel(mall);

        for inm = 1:Nm
            imm = mall(inm);

            Ylm = harmonicY( ill, imm, t3d, p3d,'type','complex','norm',true) .*(r3d.^(ill));
            Ylm = reshape(Ylm, SZ3d_c);
            
            tmp = Ylm .* dzBz; % * CHI * (-1/MUZ);
            Mcoeff = DX*DY*DZ * sum(tmp(:));

            coeffs(ifn+1, inl, inm) = Mcoeff * sqrt(4*pi/(2*ill+1));


        end
    end
    figure(10); plotsexy( (0:Nfn)/20, real(squeeze(coeffs(:,2,2))),'b-o');
    figure(30); plotsexy( (0:Nfn)/20, real(squeeze(coeffs(:,4,4))),'b-o');
    figure(32); plotsexy( (0:Nfn)/20, sqrt(2)*real(squeeze(coeffs(:,4,6))),'b-o');
    figure(50); plotsexy( (0:Nfn)/20, real(squeeze(coeffs(:,6,6))),'b-o');
    figure(52); plotsexy( (0:Nfn)/20, sqrt(2)*real(squeeze(coeffs(:,6,8))),'b-o');
    figure(54); plotsexy( (0:Nfn)/20, sqrt(2)*real(squeeze(coeffs(:,6,10))),'b-o');
    drawnow


end
save('./matfiles_new/coeffs_Mz_mur1p05_190311.mat','coeffs');

figure(21); imagesc(squeeze(dzBz(101,:,:))); axis image; colormap(cm_rb);
figure(22); imagesc(squeeze(dzBz(:,5,:))); axis image; colormap(cm_rb);


figure(121); imagesc(squeeze(Bz_rs(101,:,:))); axis image; colormap(cm_rb);
figure(122); imagesc(squeeze(Bz_rs(:,5,:))); axis image; colormap(cm_rb);
