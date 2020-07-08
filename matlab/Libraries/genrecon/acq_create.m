function [ acq_struct ] = acq_create( ro_struct, Nro, fmin, fmax, xrecon )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
    %% create sampling pattern


    acq_struct.Nro = Nro;
    acq_struct.freqs = linspace(fmin, fmax, acq_struct.Nro+1);
    acq_struct.freqs = acq_struct.freqs(1:end-1);
    acq_struct.ro_struct = ro_struct;
    acq_struct.xrecon = xrecon;
    
    %%% aliased versions of each other...
    acq_struct.fmin = fmin;
    acq_struct.fmax = fmax; 
    
    BW = fmax-fmin;
    Tsamp = 1/BW;
    
    acq_struct.tsamp = Tsamp;
    %% create recon basis

%     rec001.Nrc = 128;
%     rec001.xpos = linspace( -0.05, 0.05, rec001.Nrc );

    %%% find x-coordinates corresponding to f-samples
    %%%  ie: generate E-matrix

    % a001_E = zeros( 


    x_tmp = zeros(1,0);
    Enc_tmp = zeros(acq_struct.Nro,0);
    ixro = 1;

    for iss = 1:acq_struct.Nro

        facq = acq_struct.freqs(iss);

        fro_tmp = ro_struct.freqs - facq;
        %%% go through all positions in RO field
        %%%   if we flip sign of frequency relative to facq, the there is an
        %%%   x-point corresponding to this frequence
        for ixx = 2:numel(ro_struct.xpos)

            if sign(fro_tmp(ixx))==sign(fro_tmp(ixx-1)) && fro_tmp(ixx)~=0 && fro_tmp(ixx-1)~=0
                continue
            elseif fro_tmp(ixx-1)==0
                continue
            elseif fro_tmp(ixx)==0
                xzero = ro_struct.xpos(ixx);
                wEnc = 1/ro_struct.dfdx(ixx);

            else
                xzero = abs(fro_tmp(ixx))   / (abs(fro_tmp(ixx))+abs(fro_tmp(ixx-1))) * ro_struct.xpos(ixx-1) + ...
                        abs(fro_tmp(ixx-1)) / (abs(fro_tmp(ixx))+abs(fro_tmp(ixx-1))) * ro_struct.xpos(ixx) ;
                wEnc  = abs(fro_tmp(ixx))   / (abs(fro_tmp(ixx))+abs(fro_tmp(ixx-1))) * (ro_struct.dfdx(ixx-1)) + ...
                        abs(fro_tmp(ixx-1)) / (abs(fro_tmp(ixx))+abs(fro_tmp(ixx-1))) * (ro_struct.dfdx(ixx)) ;
            end

            x_tmp = [ x_tmp, xzero ];


            wEnc_vec = zeros( acq_struct.Nro,1);
            wEnc_vec(iss) = abs(1/wEnc);

            Enc_tmp = cat(2, Enc_tmp, wEnc_vec);

        end


    end
    acq_struct.Enc = Enc_tmp;
    acq_struct.xacq = x_tmp;

    %% create full-rank version of encoding matrix (ie no aliasing)

    iss_noalias = zeros(0,1);
    ixx_noalias = zeros(0,1);

    for iss = 1:acq_struct.Nro

        row_tmp = acq_struct.Enc(iss,:);
        iinz = find(row_tmp);
        if numel( iinz )==0 || numel( iinz )>1
            continue
        end
        iss_noalias = cat(1, iss_noalias, iss);
        ixx_noalias = cat(1, ixx_noalias, iinz);


    end

    acq_struct.Enc_noalias = acq_struct.Enc( iss_noalias, ixx_noalias);
    acq_struct.xacq_noalias = acq_struct.xacq( ixx_noalias);
    acq_struct.ind_f_noalias = iss_noalias;
    acq_struct.ind_x_noalias = ixx_noalias;

    figure; imagesc( acq_struct.Enc );
    colorbar
    axis image
    colormap gray

    figure; imagesc( acq_struct.Enc_noalias );
    colorbar
    axis image
    colormap gray

end

