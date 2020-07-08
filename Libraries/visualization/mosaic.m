function [  ] = mosaic( imgs, row_num, col_num, fig_num, title_str, disp_range, bgcolor )

% cat 2d images together as given in row_num, col_num

% convert to abs if input is complex

imag_part = imag(imgs(:));
if norm(imag_part) ~=0 
    imgs = abs(imgs);
end
    

if row_num * col_num ~= size(imgs,3)
    
    if row_num * col_num > size(imgs,3)
        % zero pad the image
        img_add = zeros(size(imgs(:,:,1)));
        imgs = cat(3, imgs, repmat(img_add,[1,1,row_num * col_num - size(imgs,3)]) ); 
    end
        
end

show_res = zeros([size(imgs,1)*row_num, size(imgs,2)*col_num]);


for r = 1:row_num
    S = imgs(:,:,col_num*(r-1)+1);
    
    for c = 2:col_num
        S = cat(2, S, imgs(:,:,col_num*(r-1)+c));
    end
    
    if r == 1
        show_res = S;
    else
        show_res = cat(1, show_res, S);
    end
end


if nargin < 6 || isempty(disp_range)
    disp_range(1) = min(show_res(:));
    disp_range(2) = max(show_res(:));
end

if nargin < 4
    figure, imagesc(show_res, disp_range), axis image off, colormap gray
else
    figure(fig_num), imagesc(show_res, disp_range), axis image off, colormap gray
end

set(gcf, 'color', 'k')
 
if nargin >= 5 && ~isempty(title_str)
    title(title_str, 'color', 'w', 'fontsize', 24)
end

drawnow


end

