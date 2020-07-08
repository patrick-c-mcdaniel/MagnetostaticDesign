function [ output ] = addBr_block(Brmat,model,Nrung, Nperlayer,domain_blk,domain_sphere2,rung_rot_angs_rad)

    Br1 = 1.32; Br2 = 1.42;  
    
    ones_ind = find(Brmat == 1);    
    twos_ind = find(Brmat == 2);
    for ii =1:Nrung/2
        rungsii_ind = [(ii-1)*Nperlayer+1:Nperlayer*ii,((ii+Nrung/2)-1)*Nperlayer+1:Nperlayer*(ii+Nrung/2)];
        ones_and_rungs_ind = intersect(ones_ind, rungsii_ind);
        twos_and_rungs_ind = intersect(twos_ind, rungsii_ind);

        %%selections groups for rungs
        selname1 = ['ones_rungs',num2str(ii),'and',num2str(ii+Nrung/2-1)];
        selname2 = ['twos_rungs',num2str(ii),'and',num2str(ii+Nrung/2-1)];

        %     model.selection.create(selname);
        model.selection(selname1).set([domain_blk(ones_and_rungs_ind)]);
        
        model.selection(selname2).set([domain_blk(twos_and_rungs_ind)]);

        %%assign to magnetic flux conservation feature
        featname1 = ['mfnc1_',num2str(ii)];
        featname2 = ['mfnc2_',num2str(ii)];

        %     model.physics('mfnc').feature.create(featname, 'MagneticFluxConservation', 3);
        model.physics('mfnc').feature(featname1).selection.named(selname1);
        model.physics('mfnc').feature(featname2).selection.named(selname2);

                
        %%set Br
        %Brnorm = 1;
        Bry1 = num2str(Br1*sin(rung_rot_angs_rad(ii)));
        Brz1 = num2str(Br1*cos(rung_rot_angs_rad(ii)));
        
        Bry2 = num2str(Br2*sin(rung_rot_angs_rad(ii)));
        Brz2 = num2str(Br2*cos(rung_rot_angs_rad(ii)));
        
        model.physics('mfnc').feature(featname1).set('Br', {'0'; Bry1; Brz1});
        model.physics('mfnc').feature(featname1).set('ConstitutiveRelationH', 'RemanentFluxDensity');
        model.physics('mfnc').feature(featname1).name([featname1 selname1]);
        
        model.physics('mfnc').feature(featname2).set('Br', {'0'; Bry2; Brz2});
        model.physics('mfnc').feature(featname2).set('ConstitutiveRelationH', 'RemanentFluxDensity');
        model.physics('mfnc').feature(featname2).name([featname2 selname2]);
    end
    
    
    %find all zeros in Brmat0
    %%overwrite default Brs with zeros in air locations
    zeros_ind  = Brmat == 0;
    selname = 'Br_zero_domains';
    zero_domains = domain_blk(zeros_ind);
    % model.selection.create(selname);
    model.selection(selname).set(zero_domains);
    
    %%assign to magnetic flux conservation feature
    featname = ['mfnc0'];
    % model.physics('mfnc').feature.create(featname, 'MagneticFluxConservation', 3);
    model.physics('mfnc').feature(featname).selection.named(selname);
    %%set Br
    model.physics('mfnc').feature(featname).set('Br', {'0'; '0'; '0'});
    model.physics('mfnc').feature(featname).set('ConstitutiveRelationH', 'RemanentFluxDensity');
    model.physics('mfnc').feature(featname).name([featname selname]);
    
    
    model.physics('mfnc').feature(featname).set('mur_mat', 'userdef');
    model.physics('mfnc').feature(featname).label('mfnc0Br_zero_domains');
    model.physics('mfnc').feature(featname).set('mur', {'1'});
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    output = 1;
%   mphsave(model, 'lastmodel');
%   
% 
% 
%     
%     
%     model.sol('sol1').runAll;
%     %
%     %
% %     model.result.numerical('av1').setResult;
% %     average_freq = mphtable(model,'tbl2');
% %     Bave = average_freq.data;
% 
%     % xyzgrid = ptgrid([-.1:.01:.1]);
%     % data = mphinterp(model,'mfnc.normB','coord',xyzgrid,'dataset','dset1');
%     datax = mpheval(model,'mfnc.Bx','selection',domain_sphere2);
%     datay = mpheval(model,'mfnc.By','selection',domain_sphere2);
%     dataz = mpheval(model,'mfnc.Bz','selection',domain_sphere2);
%     datanorm = mpheval(model,'mfnc.normB','selection',domain_sphere2);
%     [coeffs_x,polyfun,Brange] = poly_coeff_comsol(datax,0,N_order);% pause(.1);
%     [coeffs_y,polyfun,Brange] = poly_coeff_comsol(datay,0,N_order);% pause(.1);
%     [coeffs_z,polyfun,Brange] = poly_coeff_comsol(dataz,0,N_order);% pause(.1);
%     [coeffs_norm,polyfun,Brange] = poly_coeff_comsol(datanorm,0,N_order);% pause(.1);


 
