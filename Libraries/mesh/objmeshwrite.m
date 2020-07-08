function [ output ] = objmeshwrite( ver, fac, fname )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
    output = 0;
    
    header = '# Blender v2.69 (sub 0) OBJ File: '''' \n# www.blender.org';
    
    unix(['touch ' fname]);
    f = fopen(fname,'w');
    fprintf(f, header);
    
    for ivv = 1:size(ver,1)
        fprintf(f, ['\nv ' num2str(ver(ivv,1),'%06f') ' ' num2str(ver(ivv,2),'%06f') ' ' num2str(ver(ivv,3),'%06f')]);
    end
    
    fprintf(f, '\ns off');
    
    for iff = 1:size(fac,1)
        fprintf(f, ['\nf ' num2str(fac(iff,1),'%u') ' ' num2str(fac(iff,2),'%u') ' ' num2str(fac(iff,3),'%u')]);
    end  
    fclose(f);

end

