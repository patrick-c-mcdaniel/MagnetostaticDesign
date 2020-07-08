function [ cost_uni, Bmag ] = cost_Buniform_rms( Bxyz, BMEAN, Bmag )
% cost_Buniform - computes cost function, where cost is given by:
%                   range({|B|}) / mean({|B|})
% Inputs:
%   Bxyz - magnetization vector, in form required by Dipole matrix (ie.
%           Mxyz(1:3) = [Mx1, My1, Mz1]; Mxyz(4:6) = [Mx2, My2, Mz2], etc)
%   BMEAN - fixed target mean magnetic field
%   Bmag (optional) - magnitude of B at each point
%
% Outputs:
%   cost_uni - cost function = (max(Bmag(:)) - min(Bmag(:))) / mean(Bmag);
%

    if nargin<3
        Bmag = sqrt((Bxyz(1:3:end).^2 + Bxyz(2:3:end).^2 + Bxyz(3:3:end).^2));
    end
    
    cost_uni = norm( Bmag - BMEAN, 2) / numel(Bmag);
    
end

