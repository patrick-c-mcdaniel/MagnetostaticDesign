tmp = ls('./COMSOLfiles/DataFiles','-1');
tmp2 = textscan(tmp,'%s');
fn_all = tmp2{1};
Nfn = numel(fn_all);

muBm_all = zeros(Nfn,1);
delB0_all = zeros(Nfn,1);

for ifn = 1:Nfn
    fn_ii = fn_all{ifn};
    %% load comsol-derived and BS-derived 3D field maps
    %   FOV = 25cm (iso)
    %   res = 5mm (iso)
    %   matrix = 51x51x51

    %%% comsol data file

    fn_comsoldat = [ './COMSOLfiles/DataFiles/' fn_ii ];

    dat_comsol = dlmread(fn_comsoldat);

    SZCM = [43 43 43];

    Bx_cm = reshape(dat_comsol(:,4), SZCM);
    By_cm = reshape(dat_comsol(:,5), SZCM);
    Bz_cm = reshape(dat_comsol(:,6), SZCM);
    Bm_cm = sqrt( Bx_cm.^2 + By_cm.^2 + Bz_cm.^2 );

    Bxyzm_cm = cat(4, Bx_cm, By_cm, Bz_cm, Bm_cm);

    Btht_cm = atan2( sqrt(Bx_cm.^2 + By_cm.^2), Bz_cm);

%     %%% biot-savart data file
% 
%     fn_bs = '/home/patmcd/Documents/HelmetMagnet/Magnet_v2/180926_v01/BSfiles/Magnet_181007_v02_1in_dat3d_7mm.txt';
% 
%     dat_bs = dlmread(fn_bs);
% 
%     SZBS = [43 43 43];
% 
%     Bx_bs = reshape(dat_bs(:,4), SZBS);
%     By_bs = reshape(dat_bs(:,5), SZBS);
%     Bz_bs = reshape(dat_bs(:,6), SZBS);
%     Bm_bs = sqrt( Bx_bs.^2 + By_bs.^2 + Bz_bs.^2 );
% 
%     Bxyzm_bs = cat(4, Bx_bs, By_bs, Bz_bs, Bm_bs);

    %% load ROI information

    %%% load 5mm-iso ROI mask (not the one used for optimization!)
    %%%   Note: (190605) - actually this is a 7mm ROI mask (data are also
    %%%   computed on a 7mm grid).

    roi_7mm = load('../../precomp/MagROI_190528_v01.mat');

    ROImsk_7mm = roi_7mm.ROImsk_or_vol;
    ROImsk_7mm_4d = repmat(ROImsk_7mm,[ 1 1 1 4 ]);

    % ROImsk_5mm_xfl = flip(ROImsk_5mm,1);
    ROImsk_7mm_zfl = flip(ROImsk_7mm,3);

    %% plot masked comsol B0 map

    Bxyz_cm_msk = Bxyzm_cm .* ROImsk_7mm_4d;

    % figure(1); dispBmvec( Bm_cm(ROImsk_5mm_zfl), ROImsk_5mm_zfl ); %caxis([0.083 0.085]);
    load('cm_rb.mat');

    muBm = mean(Bm_cm(ROImsk_7mm));
    dispBmvec( Bm_cm(ROImsk_7mm),     ROImsk_7mm, 'mag', 'XZ', 2 ); colormap(cm_rb); caxis(muBm+[-5e-4 5e-4]);

    % dispBmvec( Btht_cm(ROImsk_5mm),     ROImsk_5mm ); title('B0 angle, comsol simulation');

    Bmean_cm = mean(Bm_cm(ROImsk_7mm));
    Bcmr = 1.5e-3 / 2;

    % dispBmvec( Bm_cm(ROImsk_5mm_zfl), ROImsk_5mm_zfl ); caxis([Bmean_cm-Bcmr, Bmean_cm+Bcmr]);
    dispBmvec( Bm_cm(ones(size(Bm_cm))>0), ones(size(Bm_cm))>0, 'mag', 'XZ', 5 ); caxis([Bmean_cm-Bcmr, Bmean_cm+Bcmr]);

    % save('./matfiles/Bm_cm_190528_v03_discreteN251.mat','Bm_cm');

    tmp = Bm_cm(ROImsk_7mm);
    delB0ROI = max(tmp)-min(tmp);
    figure(11); hist(tmp,20); grid on; grid minor;
    
    muBm_all(ifn) = muBm;
    delB0_all(ifn) = delB0ROI;
    
    disp( [ fn_ii, ' ---- B0(mean) = ', num2str(muBm), '  delB0 = ', num2str(delB0ROI) ]);
end