function [ LENTOT ] = writeBScont_xyz( cont_xyz, BSFN, Icoil, Icoil_name )
% Write Biot-Savart file for a set of wire paths
%
% Input Variables:
%
%   cont_xyz        : {Ncont}-size cell array. Each entry is a (3)x(Np) set
%                       matrix containing the set of points for one contour 
%   BSFN            : Filename (including extension) to save the BS-file
%
% Output Variables:
%
%   LENTOT          : Total combined length of all wire paths

    if nargin<3
        Icoil=1;
        Icoil_name = 'Icoil';
    end

    Nloop = numel(cont_xyz);

    unix(['touch ' BSFN]);

    fobj = fopen(BSFN,'w');

    str_head = [ 'Info {BiotSavart 4.1 data file}\n\nCurrent {\nname {Current}\nsupplies {\n	{ {' Icoil_name '} {' num2str(Icoil) '} }\n}\n}\n\n' ] ;

    fprintf(fobj, str_head);
    LENTOT = 0;

    for ilp = 1:Nloop

        tmp_str = ['Wire {\nname {Wire ' num2str(ilp) '}\ncolor 19660 45874 45874\ncurrentSupply ' Icoil_name '\nwireDiameter 0\nwinding 1\npath {\n'];
        fprintf(fobj, tmp_str);


        lptmp = cont_xyz{ilp};
        tmp_str = ['start ' num2str(lptmp(1,1), '%.5f') ' ' num2str(lptmp(2,1), '%.5f') ' ' num2str(lptmp(3,1), '%.5f') '\n'];
        fprintf(fobj, tmp_str);

        for ipts = 2:size(lptmp,2)

            tmp_str = ['goto ' num2str(lptmp(1,ipts), '%.5f') ' ' num2str(lptmp(2,ipts), '%.5f') ' ' num2str(lptmp(3,ipts), '%.5f') '\n'];
            fprintf(fobj, tmp_str);

            LENTOT = LENTOT + sqrt( (lptmp(1,ipts)-lptmp(1,ipts-1))^2 + (lptmp(2,ipts)-lptmp(2,ipts-1))^2 + (lptmp(3,ipts)-lptmp(3,ipts-1))^2 );
        end
        tmp_str = ['goto ' num2str(lptmp(1,1), '%.5f') ' ' num2str(lptmp(2,1), '%.5f') ' ' num2str(lptmp(3,1), '%.5f') '\n'];
        fprintf(fobj, tmp_str);

        LENTOT = LENTOT + sqrt( (lptmp(1,end)-lptmp(1,1))^2 + (lptmp(2,end)-lptmp(2,1))^2 + (lptmp(3,end)-lptmp(3,1))^2 );
        tmp_str = ['}\nfluxPathStep 0.1\n}\n'];
        fprintf(fobj, tmp_str);

    end

    fclose(fobj);



end