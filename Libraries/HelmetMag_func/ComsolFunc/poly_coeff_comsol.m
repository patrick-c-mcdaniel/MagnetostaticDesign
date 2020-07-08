function  [C,poly_func,Brange] = poly_coeff_comsol(data,plotswitch,N_order)

%%% defines points to generate polynomial components%%%%%% 
% delta_pos = .01;    %space between positions in field measurement (m)
% max_pos = .1;      %maximum position in X, Y, and Z (m)
% Nsamp_pos = 2*max_pos/delta_pos+1;
% [X,Y,Z] = meshgrid(-max_pos:delta_pos:max_pos, -max_pos:delta_pos:max_pos, -max_pos:delta_pos:max_pos); %3D grid of sampled field positions
% pos(:,1) = reshape(X,1, []); pos(:,2) = reshape(Y,1,[]); pos(:,3) = reshape(Z,1,[]);  %reorder sample position to array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% use points from simulated field to generate polynomial components%%%%%% 
pos(:,1) = data.p(1,:).'; pos(:,2) = data.p(2,:).'; pos(:,3) = data.p(3,:).';



field = data.d1.';

Brange = 42.58*(max(field)-min(field))*1000;

%N_order = 4; %highest polynomial order to fit

ii=1;    
 for cc= 0:(N_order)  %%% z order     
    for bb=0:(N_order)  %%% y order      
        for aa=0:(N_order) %%% x order         
            if aa+bb+cc <= N_order  %% only calculate up to N_order
                  A(:,ii) = (pos(:,1).^aa).*(pos(:,2).^bb).*(pos(:,3).^cc);   %%generates unit coefficient polynomial functions (vector form)
                  poly_func{ii} = ['X',num2str(aa),'Y',num2str(bb),'Z',num2str(cc)]; %%keeps track of polynomial functions calculated
                ii=ii+1;   
            end
        end
    end
 end




% A3D = reshape(A,Nsamp_pos,Nsamp_pos,Nsamp_pos,ii-1);  %%3D matrix version
 

% solve for coefficients from measured/simulated field
 C = A\field;
 
if plotswitch == 1

%% regrid and display 3D measured field and fit field
res = .005;
xyi = -.1:res:.1;
xyi2 = fliplr(xyi);
[qx,qy,qz] = ndgrid(xyi,xyi,xyi);
F = TriScatteredInterp(pos(:,1),pos(:,2),pos(:,3),field);
F2 = TriScatteredInterp(pos(:,1),pos(:,2),pos(:,3),A*C);

qB = 42.58*F(qx,qy,qz);
qBfit = 42.58*F2(qx,qy,qz);

% figure(210);
% subplot(1,2,1);
% image3d(qx,qy,qz,qB,1); title('simulated field');
% subplot(1,2,2);
% image3d(qx,qy,qz,qBfit,1); title('polynomial fit field');
% 
figure(111);
subplot(2,3,1);
imagesc(xyi,xyi2,squeeze(qB(21,:,:))); title('simulated field'); axis square;%caxis([3,3.108]);
xlabel('z'); ylabel('y');
subplot(2,3,2);
imagesc(xyi,xyi2,rot90(squeeze(qB(:,21,:)))); title('simulated field'); axis square;%caxis([3,3.108]);
xlabel('x'); ylabel('z');
subplot(2,3,3);
imagesc(xyi,xyi2,rot90(squeeze(qB(:,:,21)))); title('simulated field'); axis square;%caxis([3,3.108]);
xlabel('x'); ylabel('y');


subplot(2,3,4);
imagesc(xyi,xyi2,squeeze(qBfit(21,:,:))); title('polynomial fit field'); axis square;%caxis([3,3.108]);
xlabel('z'); ylabel('y');
subplot(2,3,5);
imagesc(xyi,xyi2,rot90(squeeze(qBfit(:,21,:)))); title('polynomial fit field'); axis square;% caxis([3,3.108]);
xlabel('x'); ylabel('z');
subplot(2,3,6);
imagesc(xyi,xyi2,rot90(squeeze(qBfit(:,:,21)))); title('polynomial fit field'); axis square;%caxis([3,3.108]);
 xlabel('x'); ylabel('y'); pause(0.2);
end
 
