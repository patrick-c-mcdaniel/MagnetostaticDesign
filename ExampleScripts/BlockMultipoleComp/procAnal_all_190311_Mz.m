Nfn = 20;
CHI = 0.05;
MUZ = 4*pi*1e-7;

Br = 1.42;
Mmat = Br / MUZ;

lall = 0:5;
Nl=numel(lall);

%% proc numerical

c1p05 = load('./matfiles_new/coeffs_Mz_mur1p05_190311.mat');
coeffs1p05 = c1p05.coeffs;
    figure(10); hold off; plotsexy( (0:Nfn)/20, real(squeeze(coeffs1p05(:,2,2))),'b-o');
    figure(30); hold off; plotsexy( (0:Nfn)/20, real(squeeze(coeffs1p05(:,4,4))),'b-o');
    figure(32); hold off; plotsexy( (0:Nfn)/20, sqrt(2)*real(squeeze(coeffs1p05(:,4,6))),'b-o');
    figure(50); hold off; plotsexy( (0:Nfn)/20, real(squeeze(coeffs1p05(:,6,6))),'b-o');
    figure(52); hold off; plotsexy( (0:Nfn)/20, sqrt(2)*real(squeeze(coeffs1p05(:,6,8))),'b-o');
    figure(54); hold off; plotsexy( (0:Nfn)/20, sqrt(2)*real(squeeze(coeffs1p05(:,6,10))),'b-o');
%     clear coeffs;
    %% proc numerical mur=0

% c1p00 = load('coeffs_Mz_mur1p00_190311.mat');
% coeffs1p00 = c1p00.coeffs;
% 
% load('coeffs_Mz_mur1p00_190311.mat');
%     figure(10); hold on; plotsexy( (0:Nfn)/20, real(squeeze(coeffs1p00(:,2,2))),'g-x');
%     figure(30); hold on; plotsexy( (0:Nfn)/20, real(squeeze(coeffs1p00(:,4,4))),'g-x');
%     figure(32); hold on; plotsexy( (0:Nfn)/20, sqrt(2)*real(squeeze(coeffs1p00(:,4,6))),'g-x');
%     figure(50); hold on; plotsexy( (0:Nfn)/20, real(squeeze(coeffs1p00(:,6,6))),'g-x');
%     figure(52); hold on; plotsexy( (0:Nfn)/20, sqrt(2)*real(squeeze(coeffs1p00(:,6,8))),'g-x');
%     figure(54); hold on; plotsexy( (0:Nfn)/20, sqrt(2)*real(squeeze(coeffs1p00(:,6,10))),'g-x');

%% proc numerical

coeffs_anal = zeros(Nfn+1, Nl, 2*Nl-1);

Ax = 0.0254;
Az = 0.0254;
Ay = linspace(0,0.0254,21);
 
coeffs_anal(:, 2, 2)  = sqrt(4*pi/3)  * Mmat * 1/2  * sqrt(3/pi)    * Ax*Ay*Az;
coeffs_anal(:, 4, 4)  = sqrt(4*pi/7)  * Mmat * 1/4  * sqrt(7/pi)    * ( 1/2*Ax*Ay*Az^3 - 1/4*Ax^3*Ay*Az - 1/4*Ax*Ay.^3*Az);
coeffs_anal(:, 4, 6)  = sqrt(4*pi/7)  * Mmat * 1/4  * sqrt(105/pi)  * ( 1/12*Ax^3*Ay*Az - 1/12*Ax*Ay.^3*Az );
coeffs_anal(:, 6, 6)  = sqrt(4*pi/11) * Mmat * 1/16 * sqrt(11/pi)   * ( 1/2*Ax*Ay*Az^5 - 5/6*Ax^3*Ay*Az^3 - 5/6*Ax*Ay.^3*Az^3 + 3/16*Ax^5*Ay*Az + 5/24*Ax^3*Ay.^3*Az+3/16*Ax*Ay.^5*Az);
coeffs_anal(:, 6, 8)  = sqrt(4*pi/11) * Mmat * 1/8  * sqrt(1155/pi) * ( 1/24*Ax^3*Ay*Az^3 - 1/24*Ax*Ay.^3*Az^3 - 1/80*Ax^5*Ay*Az + 1/80*Ax*Ay.^5*Az);
coeffs_anal(:, 6, 10) = sqrt(4*pi/11) * Mmat * 3/16 * sqrt(385/pi)  * ( 1/80*Ax^5*Ay*Az + 1/80*Ax*Ay.^5*Az - 1/24*Ax^3*Ay.^3*Az );

    
    figure(10); hold on; plotsexy( (0:Nfn)/20, real(squeeze(coeffs_anal(:,2,2))),'r--');
    figure(30); hold on; plotsexy( (0:Nfn)/20, real(squeeze(coeffs_anal(:,4,4))),'r--');
    figure(32); hold on; plotsexy( (0:Nfn)/20, real(squeeze(coeffs_anal(:,4,6))),'r--');
    figure(50); hold on; plotsexy( (0:Nfn)/20, real(squeeze(coeffs_anal(:,6,6))),'r--');
    figure(52); hold on; plotsexy( (0:Nfn)/20, real(squeeze(coeffs_anal(:,6,8))),'r--');
    figure(54); hold on; plotsexy( (0:Nfn)/20, real(squeeze(coeffs_anal(:,6,10))),'r--');
    drawnow
    
%% proc numerical (mur=1.05) - analytical
    figure(110); hold off; plotsexy( (0:Nfn)/20, real(squeeze(coeffs1p05(:,2,2) - coeffs_anal(:,2,2))),'b-o');
    figure(130); hold off; plotsexy( (0:Nfn)/20, real(squeeze(coeffs1p05(:,4,4) - coeffs_anal(:,4,4))),'b-o');
    figure(132); hold off; plotsexy( (0:Nfn)/20, sqrt(2)*real(squeeze(coeffs1p05(:,4,6))) - real(squeeze(coeffs_anal(:,4,6))),'b-o');
    figure(150); hold off; plotsexy( (0:Nfn)/20, real(squeeze(coeffs1p05(:,6,6) - coeffs_anal(:,6,6))),'b-o');
    figure(152); hold off; plotsexy( (0:Nfn)/20, sqrt(2)*real(squeeze(coeffs1p05(:,6,8))) - real(squeeze(coeffs_anal(:,6,8))),'b-o');
    figure(154); hold off; plotsexy( (0:Nfn)/20, sqrt(2)*real(squeeze(coeffs1p05(:,6,10))) - real(squeeze(coeffs_anal(:,6,10))),'b-o');
    
diff_1p05_L1M0 = real(squeeze(coeffs1p05(:,2,2) - coeffs_anal(:,2,2)));
diff_1p05_L3M0 = real(squeeze(coeffs1p05(:,4,4) - coeffs_anal(:,4,4)));
diff_1p05_L3M2 = sqrt(2)*real(squeeze(coeffs1p05(:,4,6))) - real(squeeze(coeffs_anal(:,4,6)));
diff_1p05_L5M0 = real(squeeze(coeffs1p05(:,6,6) - coeffs_anal(:,6,6)));
diff_1p05_L5M2 = sqrt(2)*real(squeeze(coeffs1p05(:,6,8))) - real(squeeze(coeffs_anal(:,6,8)));
diff_1p05_L5M4 = sqrt(2)*real(squeeze(coeffs1p05(:,6,10))) - real(squeeze(coeffs_anal(:,6,10)));
    
%% proc numerical (mur=1.00) - analytical
%     figure(110); hold on; plotsexy( (0:Nfn)/20, real(squeeze(coeffs1p00(:,2,2) - coeffs_anal(:,2,2))),'g-x');
%     figure(130); hold on; plotsexy( (0:Nfn)/20, real(squeeze(coeffs1p00(:,4,4) - coeffs_anal(:,4,4))),'g-x');
%     figure(132); hold on; plotsexy( (0:Nfn)/20, sqrt(2)*real(squeeze(coeffs1p00(:,4,6))) - real(squeeze(coeffs_anal(:,4,6))),'g-x');
%     figure(150); hold on; plotsexy( (0:Nfn)/20, real(squeeze(coeffs1p00(:,6,6) - coeffs_anal(:,6,6))),'g-x');
%     figure(152); hold on; plotsexy( (0:Nfn)/20, sqrt(2)*real(squeeze(coeffs1p00(:,6,8))) - real(squeeze(coeffs_anal(:,6,8))),'g-x');
%     figure(154); hold on; plotsexy( (0:Nfn)/20, sqrt(2)*real(squeeze(coeffs1p00(:,6,10))) - real(squeeze(coeffs_anal(:,6,10))),'g-x');   
% 
% diff_1p00_L1M0 = real(squeeze(coeffs1p00(:,2,2) - coeffs_anal(:,2,2)));
% diff_1p00_L3M0 = real(squeeze(coeffs1p00(:,4,4) - coeffs_anal(:,4,4)));
% diff_1p00_L3M2 = sqrt(2)*real(squeeze(coeffs1p00(:,4,6))) - real(squeeze(coeffs_anal(:,4,6)));
% diff_1p00_L5M0 = real(squeeze(coeffs1p00(:,6,6) - coeffs_anal(:,6,6)));
% diff_1p00_L5M2 = sqrt(2)*real(squeeze(coeffs1p00(:,6,8))) - real(squeeze(coeffs_anal(:,6,8)));
% diff_1p00_L5M4 = sqrt(2)*real(squeeze(coeffs1p00(:,6,10))) - real(squethe eze(coeffs_anal(:,6,10)));
    
%% poly-fit differences

pf_1p05_L1M0 = polyfit( Ay, diff_1p05_L1M0', 5);
pf_1p05_L3M0 = polyfit( Ay, diff_1p05_L3M0', 5);
pf_1p05_L3M2 = polyfit( Ay, diff_1p05_L3M2', 5);
pf_1p05_L5M0 = polyfit( Ay, diff_1p05_L5M0', 5);
pf_1p05_L5M2 = polyfit( Ay, diff_1p05_L5M2', 5);
pf_1p05_L5M4 = polyfit( Ay, diff_1p05_L5M4', 5);

diff_pf_1p05_L1M0 = polyval(pf_1p05_L1M0, Ay);
diff_pf_1p05_L3M0 = polyval(pf_1p05_L3M0, Ay);
diff_pf_1p05_L3M2 = polyval(pf_1p05_L3M2, Ay);
diff_pf_1p05_L5M0 = polyval(pf_1p05_L5M0, Ay);
diff_pf_1p05_L5M2 = polyval(pf_1p05_L5M2, Ay);
diff_pf_1p05_L5M4 = polyval(pf_1p05_L5M4, Ay);

    figure(110); hold on; plotsexy( (0:Nfn)/20, diff_pf_1p05_L1M0,'r--');
    figure(130); hold on; plotsexy( (0:Nfn)/20, diff_pf_1p05_L3M0,'r--');
    figure(132); hold on; plotsexy( (0:Nfn)/20, diff_pf_1p05_L3M2,'r--');
    figure(150); hold on; plotsexy( (0:Nfn)/20, diff_pf_1p05_L5M0,'r--');
    figure(152); hold on; plotsexy( (0:Nfn)/20, diff_pf_1p05_L5M2,'r--');
    figure(154); hold on; plotsexy( (0:Nfn)/20, diff_pf_1p05_L5M4,'r--');  
    
    save('matfiles_new/pf_1p05_L1M0.mat','pf_1p05_L1M0');
    save('matfiles_new/pf_1p05_L3M0.mat','pf_1p05_L3M0');
    save('matfiles_new/pf_1p05_L3M2.mat','pf_1p05_L3M2');
    save('matfiles_new/pf_1p05_L5M0.mat','pf_1p05_L5M0');
    save('matfiles_new/pf_1p05_L5M2.mat','pf_1p05_L5M2');
    save('matfiles_new/pf_1p05_L5M4.mat','pf_1p05_L5M4');

% save('coeffs_Mz_mur1p05_190311.mat','coeffs');
% 
% figure(21); imagesc(squeeze(dzBz(101,:,:))); axis image; colormap(cm_rb);
% figure(22); imagesc(squeeze(dzBz(:,5,:))); axis image; colormap(cm_rb);
% 
% 
% figure(121); imagesc(squeeze(Bz_rs(101,:,:))); axis image; colormap(cm_rb);
% figure(122); imagesc(squeeze(Bz_rs(:,5,:))); axis image; colormap(cm_rb);
