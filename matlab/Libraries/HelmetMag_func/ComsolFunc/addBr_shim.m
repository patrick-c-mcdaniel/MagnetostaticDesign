function [ output ] = addBr_shim(model, domain_shim, Brxyz)

    Nshim = numel(domain_shim);
    
    for ishm =1:Nshim
        disp(['adding shim Br condition, block: ' num2str(ishm) ' / ' num2str(Nshim)]);
        selname = ['shim_' num2str(ishm,'%04u') '_Brcond'];
%         if ishm~=0
%             
%             model.selection.create(selname);
%         end
%         tic
%         model.selection(selname).set(domain_shim(ishm));
       
    tic
        %%assign to magnetic flux conservation feature
        
            featname = ['mfnc_' num2str(ishm,'%04u')];
%         if ishm~=1
            model.physics('mfnc').feature.create(featname, 'MagneticFluxConservation', 3);
            model.physics('mfnc').feature(featname).selection.set(domain_shim(ishm));
%         end
        toc
                

        Brx = num2str(Brxyz(ishm,1));
        Bry = num2str(Brxyz(ishm,2));
        Brz = num2str(Brxyz(ishm,3));
        
        model.physics('mfnc').feature(featname).set('Br', {Brx; Bry; Brz});
        model.physics('mfnc').feature(featname).set('ConstitutiveRelationH', 'RemanentFluxDensity');
        model.physics('mfnc').feature(featname).name([featname selname]);
        toc
    end
    
    

    
    output = 1;
end

 
