function [ Call ] = nlcon_Mmax( Mxyz, Mmax )
% nlcon_Mmax - calculates whether magnetization vector meets nonlinear max
%              magnetization constraint or not
% Inputs:
%   Mxyz - magnetization vector, in form required by Dipole matrix (ie.
%           Mxyz(1:3) = [Mx1, My1, Mz1]; Mxyz(4:6) = [Mx2, My2, Mz2], etc)
%   Mmax - maximum allowable magnetization magnitude, in SI units of
%           magnetic dipole moment (equal to: Br[T] * Vol[m^3] / MUZ )
%
% Outputs:
%   Call - binary vector specifying whether each individual magnetization
%           meets the constraint or not
%             1==Meets constraint (|M|<=Mmax)
%             0==Does not meet constraint (|M|>Mmax)

%     Callb = (Mxyz(1:3:end).^2 + Mxyz(2:3:end).^2 + Mxyz(3:3:end).^2) < Mmax^2;
%     Canyb = sum(Callb(:))>0;
%     Cany = Canyb==0;
%     Call = Callb==0;

    Call = (Mxyz(1:3:end).^2 + Mxyz(2:3:end).^2 + Mxyz(3:3:end).^2) - Mmax^2;

end

