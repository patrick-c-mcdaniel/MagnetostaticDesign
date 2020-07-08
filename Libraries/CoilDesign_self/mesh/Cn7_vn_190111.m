function [Cx,Cy,Cz] = Cn7_vn_190111(node,triangle,basis,r,ver_vn,cell_vn)

% This file aim at calculating the so called "Cn matrix"
%
%
% node : a matrix with the 3d position of each node (in meter)
% triangle : a matrix linking 3 node together to form a triangle
% basis : structure with the basis vector and the related informations
% r : the target point

%ChangeLog
% v5: adapted to the new structure of node and triangle
% v7: clean the code

tic;


[u,v,ck] = triGaussPoints(2);
for i=1:size(u,1)
    w(i) = 1-u(i)-v(i);
end


mu_0 = 4*pi*10^-7;
coef = mu_0/(4*pi);
nbrIntegrationPoints = size(ck,1);

Cx = zeros(size(r,1),size(node,2));
Cy = zeros(size(r,1),size(node,2));
Cz = zeros(size(r,1),size(node,2));

%activate the parallel function
matlabVersion = version;
matlabVersion = str2num(matlabVersion(1:3));
if matlabVersion < 8.2
    [TF,~] = license('checkout', 'Distrib_Computing_Toolbox');
    if TF
        schd = findResource('scheduler', 'configuration', 'local');
        numWorkers = schd.ClusterSize;
    end

    if matlabpool('size') == 0  && TF && numWorkers >1
        % checking to see if the pool is already open and of we have the licence
        % and at least 2 cores
        matlabpool open
    end
elseif matlabVersion >= 8.2 
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
	if isempty(poolobj)
		parpool;
    end
end

% fprintf('Loop n=%5.0i',0);
for n=1:size(node,2) % For a given node
    tempCx = zeros(1,size(r,1));
    tempCy = zeros(1,size(r,1));
    tempCz = zeros(1,size(r,1));

    ptstmp = ver_vn( cell_vn{n},:);
    xtmp = [ptstmp(:,1); ptstmp(1,1)];
    ytmp = [ptstmp(:,2); ptstmp(1,2)];
    ztmp = [ptstmp(:,3); ptstmp(1,3)];
    %%% see which way the cell turns...

    pt1 = [ xtmp(1), ytmp(1), ztmp(1) ];
    pt2 = [ xtmp(2), ytmp(2), ztmp(2) ];
    pt3 = [ xtmp(3), ytmp(3), ztmp(3) ];

    v21 = pt2-pt1;
    v32 = pt3-pt2;

    vcross = cross( v21', v32' );
    if vcross(1)>0
        disp(['n = ' num2str(n) ' , x = (+)']);
    else
        disp(['n = ' num2str(n) ' , x = (-)']);
        xtmp = flip(xtmp);
        ytmp = flip(ytmp);
        ztmp = flip(ztmp);  
    end
    
    for k=1:size(r,1); % And for a given target point
        bigSumCx = 0;
        bigSumCy = 0;
        bigSumCz = 0;
        

        
        for iseg = 1:(numel(xtmp)-1)
            r_x = r(k,1); %X coordinate of the provided point
            r_y = r(k,2); %Y coordinate of the provided point
            r_z = r(k,3); %Z coordinate of the provided point
            
            vni_x = xtmp(iseg+1)-xtmp(iseg);
            vni_y = ytmp(iseg+1)-ytmp(iseg);
            vni_z = ztmp(iseg+1)-ztmp(iseg);
            
            
            rpk_x = (xtmp(iseg+1)+xtmp(iseg))/2;
            rpk_y = (ytmp(iseg+1)+ytmp(iseg))/2;
            rpk_z = (ztmp(iseg+1)+ztmp(iseg))/2;
            for l=1 %:nbrIntegrationPoints

                normRRp = sqrt((r_x-rpk_x(l))^2+(r_y-rpk_y(l))^2+(r_z-rpk_z(l))^2);
                
                fCx(l) = ck(l)*((-vni_z*(r_y-rpk_y(l))+(vni_y*(r_z-rpk_z(l))))/(normRRp^3));
                fCy(l) = ck(l)*((-vni_x*(r_z-rpk_z(l))+(vni_z*(r_x-rpk_x(l))))/(normRRp^3));
                fCz(l) = ck(l)*((-vni_y*(r_x-rpk_x(l))+(vni_x*(r_y-rpk_y(l))))/(normRRp^3));
            end
            
            sum1Cx = sum(fCx);
            sum1Cy = sum(fCy);
            sum1Cz = sum(fCz);
            
            bigSumCx = bigSumCx + sum1Cx;
            bigSumCy = bigSumCy + sum1Cy;
            bigSumCz = bigSumCz + sum1Cz;

        end
        tempCx(k) = coef*bigSumCx;
        tempCy(k) = coef*bigSumCy;
        tempCz(k) = coef*bigSumCz;
        %C(k,n) = coef*bigSum;
    end
    Cx(:,n) = tempCx; %*node(n).nbrTriangle;
    Cy(:,n) = tempCy; %*node(n).nbrTriangle;
    Cz(:,n) = tempCz; %*node(n).nbrTriangle;
%     fprintf(1,'\b\b\b\b\b%5.0i',n);
end
%fprintf(1,' - Done \n');
fprintf(' - Done in %5.0f sec.\n',toc);