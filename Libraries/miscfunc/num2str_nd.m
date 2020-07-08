function s = num2str_nd(varargin)
    tmp = num2str(varargin{:});
    tmp2 = tmp;
    for ichar = 1:numel(tmp)
       
        if strcmp(tmp(ichar),'.')
            tmp2(ichar) = 'p';
        else
            tmp2(ichar) = tmp(ichar);
        end
        
    end
    s = tmp2;

end