function [ s_t, tall, s_f ] = acq_signal( a_struct, Tapod )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
    
    if nargin<2
        Tapod = Inf;
    end
   
    
    
    %% acq. data

    phant = a_struct.phant_struct;
    xphant = phant.xpts;
    
    fmax = a_struct.fmax;
    fmin = a_struct.fmin;
    
    BW = fmax-fmin;
    Tsamp = 1/BW;
    Fsamp = (a_struct.freqs(2) - a_struct.freqs(1));
    
    tall = Tsamp * ((-a_struct.Nro/2+1/2):(a_struct.Nro/2-1/2));
    
    
    delXphant = xphant(2) - xphant(1);
    
    fphant = interp1( a_struct.ro_struct.xpos, a_struct.ro_struct.freqs, xphant, 'cubic' );
    
    f_2d = repmat( fphant(:), [1, a_struct.Nro] );
    M_2d = repmat( phant.Mag(:), [1, a_struct.Nro] );
    
    t_2d = repmat( tall(:)', [phant.Npts, 1] );
    
    GYRO = 42.576e6;
    
    s_t_2d = M_2d .* exp(1j*2*pi*f_2d.*t_2d) ; 
    s_t = Fsamp*sum( delXphant * s_t_2d, 1);
    
    %% apod. data
    
    W_apod = exp( -1* ( abs(tall)/Tapod));
%    W_apod = hamming( numel(tall));
%      s_t = W_apod(:) .* s_t(:);
    
    s_f = fft( s_t(:) );
%     
%     HzPx = (Fsamp) / a_struct.ro_struct.delx;
%     
    

end

