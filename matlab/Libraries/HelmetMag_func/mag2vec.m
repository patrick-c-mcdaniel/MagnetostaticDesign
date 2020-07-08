function [ mvec ] = mag2vec( mx, my, mz )
% mag2vec - turns lists of x-, y-, and z-components of magnetization into a
% single vector, properly arranged for forward Dipole matrix calculation

%     tmp = cat(2, mx(:), my(:), mz(:));
%     tmp = tmp';
%     mvec = tmp(:);
    mvec = zeros(numel(mx)+numel(my)+numel(mz), 1);
    mvec(1:3:end) = mx;
    mvec(2:3:end) = my;
    mvec(3:3:end) = mz;
    

end

