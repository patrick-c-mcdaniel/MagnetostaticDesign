function [ rgball ] = dat2colors( datvec, cmap ) 
    %% function [ rbgall ] = dat2colors( datvec, cmap ) 
    %     Create list of rgb colors given data vector and discrete colormap
    %       Automatically scales the colormap to match the datarange

    dat_min = min(datvec);
    dat_max = max(datvec);
    Ncm = size(cmap,1);
    
    DELdat = (dat_max-dat_min)/(Ncm-1);
    
    ind_all = round((datvec - dat_min)/DELdat + 1);
    
    rgball = cmap( ind_all(:),:);
    
end
