function [phn] = imagesc_bbg(varargin)
    

    maskROI = varargin{1};
    
    image( zeros([size(varargin{2}), 3])); hold on;
    phn = imagesc(varargin{2:end}, 'AlphaData',double(maskROI));
    
    

end