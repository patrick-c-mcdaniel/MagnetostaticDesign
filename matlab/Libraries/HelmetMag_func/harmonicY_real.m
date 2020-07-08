function [ sh_r ] = harmonicY_real( varargin )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


%%% ill, imm, tt3, pp3,'type','complex','norm',true) .*(rr3.^(ill));
    m = varargin{2};
    varargin_neg = varargin;
    varargin_neg{2} = -1*m;
    
    if m==0
        sh_r = harmonicY( varargin{:} );
    elseif m<0
        sh_r = 1j/sqrt(2) * (harmonicY( varargin{:} )     - (-1)^m*harmonicY( varargin_neg{:} ));
    elseif m>0
        sh_r =  1/sqrt(2) * (harmonicY( varargin_neg{:} ) + (-1)^m*harmonicY( varargin{:} ));
    end
end

