function out = model
%
% Block_mur_190225_matlab.m
%
% Model exported on Feb 26 2019, 16:19 by COMSOL 5.1.0.180.

import com.comsol.model.*
import com.comsol.model.util.*

%% loop through block sizes (50 values from 0.02in to 1.00in)
delY = 0.0254/20;
y_all = delY:delY:(20*delY);

for iyy = 1:numel(y_all)
        
    disp(['iyy = ' num2str(iyy) ' / ' num2str(20) ]);
    szy = y_all(iyy);
    szy_str = num2str(szy,10);
    
    mesh_min = szy/10;
    mesh_max = szy/10*2;
    
    mesh_min_str = num2str(mesh_min,10);
    mesh_max_str = num2str(mesh_max,10);
    
    %%% y sampling points
    yexp_del = delY/10;
    yexp_min = -1*( iyy*5-1 )*yexp_del - yexp_del/2;
    yexp_max =  1*( iyy*5-1 )*yexp_del + yexp_del/2;
    
    yexp_del_str = num2str(yexp_del,10);
    yexp_min_str = num2str(yexp_min,10);
    yexp_max_str = num2str(yexp_max,10);
    
    %%% x sampling points
    xexp_min = -1*( 99 )*yexp_del - yexp_del/2;
    xexp_max =  1*( 99 )*yexp_del + yexp_del/2;
    
    xexp_del_str = num2str(yexp_del,10);
    xexp_min_str = num2str(xexp_min,10);
    xexp_max_str = num2str(xexp_max,10);
    
    %%% z sampling points
    zexp_min = -1*( 99 )*yexp_del - 7*yexp_del/2;
    zexp_max =  1*( 99 )*yexp_del + 7*yexp_del/2;
   
    zexp_del_str = num2str(yexp_del,10);
    zexp_min_str = num2str(zexp_min,10);
    zexp_max_str = num2str(zexp_max,10);
    
    fn_datexport = ['Data/Block_y' num2str(iyy,'%02u') '_Bxyz_190311_Mz.txt'];

    model = ModelUtil.create('Model');

    model.modelPath('./Comsol');

    model.comments(['Untitled\n\n']);

    model.modelNode.create('comp1');

    model.geom.create('geom1', 3);

    model.mesh.create('mesh1', 'geom1');

    model.physics.create('mfnc', 'MagnetostaticsNoCurrents', 'geom1');

    model.geom('geom1').create('blk1', 'Block');
    model.geom('geom1').feature('blk1').label('Block_n01');
    model.geom('geom1').feature('blk1').set('size', {'0.0254' '0.0254' '1'});
    model.geom('geom1').feature('blk1').setIndex('size', '0.0254', 2);
    model.geom('geom1').feature('blk1').setIndex('size', szy_str, 1);
    model.geom('geom1').feature('blk1').set('base', 'center');
    model.geom('geom1').run('blk1');
    model.geom('geom1').run('blk1');
    model.geom('geom1').create('sph1', 'Sphere');
    model.geom('geom1').feature('sph1').set('r', '0.5');
    model.geom('geom1').run('sph1');
    model.geom('geom1').run;

    model.material.create('mat1', 'Common', 'comp1');
    model.material('mat1').label('NdFeB');
    model.material('mat1').set('family', 'plastic');
    model.material('mat1').propertyGroup('def').set('relpermeability', '1.05');
    model.material('mat1').propertyGroup.create('BHCurve', 'BH curve');
    model.material('mat1').propertyGroup('BHCurve').set('normB', '');
    model.material('mat1').propertyGroup('BHCurve').set('normH', 'sqrt(H1^2+H2^2+H3^2)');
    model.material('mat1').propertyGroup('BHCurve').addInput('magneticfield');
    model.material('mat1').set('family', 'plastic');
    model.material.create('mat2', 'Common', 'comp1');
    model.material('mat2').label('Air');
    model.material('mat2').set('family', 'air');
    model.material('mat2').propertyGroup('def').set('relpermeability', '1');
    model.material('mat2').propertyGroup('def').set('relpermittivity', '1');
    model.material('mat2').propertyGroup('def').set('dynamicviscosity', 'eta(T[1/K])[Pa*s]');
    model.material('mat2').propertyGroup('def').set('ratioofspecificheat', '1.4');
    model.material('mat2').propertyGroup('def').set('electricconductivity', '0[S/m]');
    model.material('mat2').propertyGroup('def').set('heatcapacity', 'Cp(T[1/K])[J/(kg*K)]');
    model.material('mat2').propertyGroup('def').set('density', 'rho(pA[1/Pa],T[1/K])[kg/m^3]');
    model.material('mat2').propertyGroup('def').set('thermalconductivity', 'k(T[1/K])[W/(m*K)]');
    model.material('mat2').propertyGroup('def').set('soundspeed', 'cs(T[1/K])[m/s]');
    model.material('mat2').propertyGroup('def').func.create('eta', 'Piecewise');
    model.material('mat2').propertyGroup('def').func('eta').set('funcname', 'eta');
    model.material('mat2').propertyGroup('def').func('eta').set('arg', 'T');
    model.material('mat2').propertyGroup('def').func('eta').set('extrap', 'constant');
    model.material('mat2').propertyGroup('def').func('eta').set('pieces', {'200.0' '1600.0' '-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4'});
    model.material('mat2').propertyGroup('def').func.create('Cp', 'Piecewise');
    model.material('mat2').propertyGroup('def').func('Cp').set('funcname', 'Cp');
    model.material('mat2').propertyGroup('def').func('Cp').set('arg', 'T');
    model.material('mat2').propertyGroup('def').func('Cp').set('extrap', 'constant');
    model.material('mat2').propertyGroup('def').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
    model.material('mat2').propertyGroup('def').func.create('rho', 'Analytic');
    model.material('mat2').propertyGroup('def').func('rho').set('funcname', 'rho');
    model.material('mat2').propertyGroup('def').func('rho').set('args', {'pA' 'T'});
    model.material('mat2').propertyGroup('def').func('rho').set('expr', 'pA*0.02897/8.314/T');
    model.material('mat2').propertyGroup('def').func('rho').set('dermethod', 'manual');
    model.material('mat2').propertyGroup('def').func('rho').set('argders', {'pA' 'd(pA*0.02897/8.314/T,pA)'; 'T' 'd(pA*0.02897/8.314/T,T)'});
    model.material('mat2').propertyGroup('def').func.create('k', 'Piecewise');
    model.material('mat2').propertyGroup('def').func('k').set('funcname', 'k');
    model.material('mat2').propertyGroup('def').func('k').set('arg', 'T');
    model.material('mat2').propertyGroup('def').func('k').set('extrap', 'constant');
    model.material('mat2').propertyGroup('def').func('k').set('pieces', {'200.0' '1600.0' '-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4'});
    model.material('mat2').propertyGroup('def').func.create('cs', 'Analytic');
    model.material('mat2').propertyGroup('def').func('cs').set('funcname', 'cs');
    model.material('mat2').propertyGroup('def').func('cs').set('args', {'T'});
    model.material('mat2').propertyGroup('def').func('cs').set('expr', 'sqrt(1.4*287*T)');
    model.material('mat2').propertyGroup('def').func('cs').set('dermethod', 'manual');
    model.material('mat2').propertyGroup('def').func('cs').set('argders', {'T' 'd(sqrt(1.4*287*T),T)'});
    model.material('mat2').propertyGroup('def').addInput('temperature');
    model.material('mat2').propertyGroup('def').addInput('pressure');
    model.material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
    model.material('mat2').propertyGroup('RefractiveIndex').set('n', '1');
    model.material('mat2').set('family', 'air');
    model.material('mat1').selection.set([2]);
    model.material('mat2').selection.set([1]);

    model.physics('mfnc').feature.create('mfc2', 'MagneticFluxConservation', 3);
    model.physics('mfnc').feature('mfc2').selection.set([2]);
    model.physics('mfnc').feature('mfc2').set('ConstitutiveRelationH', 'RemanentFluxDensity');
    model.physics('mfnc').feature('mfc2').set('Br', {'0' '0' '1.42'});

    model.mesh('mesh1').automatic(false);
    model.mesh('mesh1').feature('size').set('hauto', '3');
    model.mesh('mesh1').feature('ftet1').create('size1', 'Size');
    model.mesh('mesh1').feature('ftet1').feature('size1').set('custom', 'on');
    model.mesh('mesh1').feature('ftet1').feature('size1').set('hmaxactive', 'on');
    model.mesh('mesh1').feature('ftet1').feature('size1').set('hmax', mesh_max_str);
    model.mesh('mesh1').feature('ftet1').selection.geom('geom1', 3);
    model.mesh('mesh1').feature('ftet1').selection.set([2]);
    model.mesh('mesh1').create('size1', 'Size');
    model.mesh('mesh1').feature('size1').set('hauto', '3');
    model.mesh('mesh1').feature('size1').set('custom', 'on');
    model.mesh('mesh1').feature('size1').set('hminactive', 'on');
    model.mesh('mesh1').feature('size1').set('hmin', mesh_min_str);
%     model.mesh('mesh1').run;
    model.mesh('mesh1').feature.remove('size1');
    model.mesh('mesh1').create('ftet2', 'FreeTet');
    model.mesh('mesh1').feature('ftet2').create('size1', 'Size');
    model.mesh('mesh1').feature('ftet2').feature('size1').set('custom', 'off');
    model.mesh('mesh1').feature('ftet2').feature('size1').set('hauto', '3');
    model.mesh('mesh1').feature('ftet2').feature('size1').set('custom', 'on');
    model.mesh('mesh1').feature('ftet2').feature('size1').set('hminactive', 'on');
    model.mesh('mesh1').feature('ftet2').feature('size1').set('hmin', mesh_min_str);
    
    model.mesh('mesh1').feature('ftet1').feature('size1').set('hmaxactive', 'on');
    model.mesh('mesh1').feature('ftet1').feature('size1').set('hminactive', 'on');
    model.mesh('mesh1').feature('ftet1').feature('size1').set('hmin', '0.0001');
    model.mesh('mesh1').feature('ftet1').feature('size1').set('hmax', '0.0002');
    model.mesh('mesh1').feature('ftet1').feature('size1').set('hnarrowactive', 'off');
    model.mesh('mesh1').feature('ftet1').feature('size1').set('hminactive', 'off');
    model.mesh('mesh1').feature('ftet1').feature('size1').set('hmax', '0.0005');
    model.mesh('mesh1').feature('ftet1').feature('size1').set('hnarrowactive', 'on');
    model.mesh('mesh1').feature('ftet1').feature('size1').set('hnarrow', '0.25');
    % model.mesh('mesh1').run('ftet1');
    model.mesh('mesh1').feature('ftet1').feature('size1').set('hmax', '0.0005');
    model.mesh('mesh1').run('ftet1');
    model.mesh('mesh1').feature('ftet2').feature('size1').set('hminactive', 'off');
    model.mesh('mesh1').feature('ftet2').feature('size1').set('hgradactive', 'on');
    model.mesh('mesh1').feature('ftet2').feature('size1').set('hgrad', '1.6');
    % model.mesh('mesh1').run('ftet2');
    model.mesh('mesh1').run;

    model.study.create('std1');
    model.study('std1').create('stat', 'Stationary');
    model.study('std1').feature('stat').activate('mfnc', true);

    model.sol.create('sol1');
    model.sol('sol1').study('std1');

    model.study('std1').feature('stat').set('notlistsolnum', 1);
    model.study('std1').feature('stat').set('notsolnum', '1');
    model.study('std1').feature('stat').set('listsolnum', 1);
    model.study('std1').feature('stat').set('solnum', '1');

    model.sol('sol1').create('st1', 'StudyStep');
    model.sol('sol1').feature('st1').set('study', 'std1');
    model.sol('sol1').feature('st1').set('studystep', 'stat');
    model.sol('sol1').create('v1', 'Variables');
    model.sol('sol1').feature('v1').set('control', 'stat');
    model.sol('sol1').create('s1', 'Stationary');
    model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
    model.sol('sol1').feature('s1').create('i1', 'Iterative');
    model.sol('sol1').feature('s1').feature('i1').set('linsolver', 'cg');
    model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('prefun', 'amg');
    model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'i1');
    model.sol('sol1').feature('s1').feature.remove('fcDef');
    model.sol('sol1').attach('std1');

    model.result.create('pg1', 'PlotGroup3D');
    model.result('pg1').label('Magnetic Flux Density Norm (mfnc)');
    model.result('pg1').set('oldanalysistype', 'noneavailable');
    model.result('pg1').set('data', 'dset1');
    model.result('pg1').feature.create('mslc1', 'Multislice');
    model.result('pg1').feature('mslc1').set('oldanalysistype', 'noneavailable');
    model.result('pg1').feature('mslc1').set('expr', 'mfnc.normB');
    model.result('pg1').feature('mslc1').set('data', 'parent');

    model.sol('sol1').runAll;

    model.result('pg1').run;


    model.label('Block_mur_190225.mph');

    model.result('pg1').run;
    model.result('pg1').feature('mslc1').set('expr', 'mfnc.Mz');
    model.result('pg1').run;
    model.result('pg1').feature('mslc1').set('rangecoloractive', 'off');
    model.result('pg1').feature('mslc1').set('expr', 'd(mfnc.Mz,z)+d(mfnc.My,y)+d(mfnc.Mx,x)');
    model.result('pg1').run;
    model.result('pg1').feature('mslc1').set('expr', 'mfnc.Mz');
    model.result('pg1').run;
    model.result.export.create('data1', 'Data');
    model.result.export('data1').setIndex('expr', 'mfnc.Mz', 0);
    model.result.export('data1').set('location', 'grid');
    model.result.export('data1').set('gridx3', ['range(' xexp_min_str ',' xexp_del_str ',' xexp_max_str ')']);
    model.result.export('data1').set('gridy3', ['range(' yexp_min_str ',' yexp_del_str ',' yexp_max_str ')']);
    model.result.export('data1').set('gridz3', ['range(' zexp_min_str ',' zexp_del_str ',' zexp_max_str ')']);
    model.result.export('data1').set('filename', fn_datexport);
    model.result.export('data1').run;
    model.result('pg1').run;
end
out = model;
