function [ varout ] = patchsexy( varargin )

    varout = patch(varargin{:});
    
    set( varout, 'FaceColor',[1 0.8 0.9]);
    axis equal

end