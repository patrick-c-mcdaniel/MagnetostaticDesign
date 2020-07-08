function [ vals, inds ] = RoundToSet( invals, disc_set )
%% RoundToSet - round real-valued numbers to discrete set of numbers
%    Inputs:
%       inval   : input set of values to be round (size: Nin x 1)
%       disc-set: discrete set of values to which to round input (size: Ndisc x 1)
%
%   Outputs:
%       vals    : rounded version of invals
%       inds    : indices into disc_set for vals
    sz_in = size(invals);
    
    invals = invals(:);
    disc_set = disc_set(:);
    
    invals2d = repmat( invals, [ 1 numel(disc_set) ]);
    disc2d   = repmat( disc_set', [ numel(invals) 1 ]);
    
    [ vmin, inds ] = min( abs( invals2d - disc2d ), [], 2 );
    vals = disc_set(inds);
    
    inds = reshape(inds, sz_in);
    vals = reshape(vals, sz_in);




end