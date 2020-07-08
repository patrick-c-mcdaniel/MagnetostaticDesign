function [  ] = mosaic( imgs, row_num, col_num, fig_num, title_str, disp_range, bgcolor, bgmask )

if nargin < 8
    bgcolor = [ 0 0 0 ];
    bgmask =  zeros(size(imgs));
end

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
        
        bg_add  = zeros(size(bgmask(:,:,1)));
        bgmask  = cat(3, bgmask, repmat(bg_add,[1, 1, row_num*col_num-size(bgmask,3)]) );
    end
        
end

show_res = zeros([size(imgs,1)*row_num, size(imgs,2)*col_num]);


for r = 1:row_num
    S = imgs(:,:,col_num*(r-1)+1);
    Sb = bgmask(:,:,col_num*(r-1)+1);
    
    for c = 2:col_num
        S = cat(2, S, imgs(:,:,col_num*(r-1)+c));
        Sb = cat(2, Sb, bgmask(:,:,col_num*(r-1)+c));
    end
    
    if r == 1
        show_res = S;
        show_bg  = Sb;
    else
        show_res = cat(1, show_res, S);  
        show_bg  = cat(1, show_bg,  Sb);
    end
end


if nargin < 6 || isempty(disp_range)
    disp_range(1) = min(show_res(:));
    disp_range(2) = max(show_res(:));
end

if nargin < 4
    figure, imagesc(show_res, disp_range), axis image off, colormap gray; hold on;
    image('CData', cat(3, bgcolor(1)*ones(size(show_bg)), bgcolor(2)*ones(size(show_bg)), bgcolor(3)*ones(size(show_bg))), 'AlphaData',show_bg); hold off;
else
    figure(fig_num), imagesc(show_res, disp_range), axis image off, colormap gray
    image('CData', cat(3, bgcolor(1)*ones(size(show_bg)), bgcolor(2)*ones(size(show_bg)), bgcolor(3)*ones(size(show_bg))), 'AlphaData',show_bg); hold off;
end

set(gcf, 'color', 'w')
 
if nargin >= 5 && ~isempty(title_str)
    title(title_str, 'color', 'w', 'fontsize', 24)
end

drawnow


end

