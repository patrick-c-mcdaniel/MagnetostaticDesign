function [ cost_uni, Bmag ] = cost_Buniform_nomean_Btarg( Bt, Bmag )
% cost_Buniform - computes cost function, where cost is given by:
%                   range({|B|}) / mean({|B|})
% Inputs:
%   Bxyz - magnetization vector, in form required by Dipole matrix (ie.
%           Mxyz(1:3) = [Mx1, My1, Mz1]; Mxyz(4:6) = [Mx2, My2, Mz2], etc)
%   Bmag (optional) - magnitude of B at each point
%
% Outputs:
%   cost_uni - cost function = (max(Bmag(:)) - min(Bmag(:))) / mean(Bmag);
%

    if nargin<2
        Bmag = Bt;
    end
    
    cost_uni = (max(Bt(:)) - min(Bt(:)));
    
end

