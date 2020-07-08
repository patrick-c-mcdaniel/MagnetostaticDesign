function output = dispimgs(imarray,cmap, orient,mincol, maxcol, voxdim, axstr)
if nargin < 7
    axstr = 'image';
end
if nargin < 6 || isequal(voxdim,[]) || isequal(voxdim,[])
    voxdim = [1 1 1];
end
if nargin < 5 || isequal(mincol,[]) || isequal(maxcol,[])
    mincol = -1;
    maxcol = -1;
end
if nargin < 3
    orient = 'axial';
end
if nargin < 2
    cmap = 'gray';
end


hfig = figure;
clf(hfig)
imarrayc = cell(size(imarray,3),1);
imarrayd = zeros(size(imarray));
if strcmp(class(imarray), 'double') || strcmp(class(imarray),'logical') 
    imarrayd = imarray;
    for i = 1:size(imarray,3)
        imarrayc{i} = imarray(:,:,i);
    end
elseif strcmp(class(imarray),'cell')
    imarrayc = imarray;
    for i = 1:numel(imarray)
        imarrayd(:,:,i) = imarray{i};
    end
end
switch orient
    case 'axial'
        smax = size(imarrayd,3);
        pixdim = [voxdim(1) voxdim(2) 1];
    case 'sag'
        smax = size(imarrayd,1);
        for i = 1:size(imarrayd,1)
            imarrayc{i} = squeeze(imarrayd(i,:,:));
        end
        pixdim = [voxdim(2) voxdim(3) 1];
    case 'cor'
        smax = size(imarrayd,2);
        for i = 1:size(imarrayd,2)
            imarrayc{i} = squeeze(imarrayd(:,i,:));
        end
        pixdim = [voxdim(1) voxdim(3) 1];
    otherwise
        smax = size(imarrayd,3);
        pixdim = [voxdim(1) voxdim(2) 1];
        print('invalid orientation; defaulting to axial');
end

if mincol == -1
    %tmp = imarrayc{:};
    mincol = min(imarrayd(:));
    maxcol = max(imarrayd(:));
end




imagesc(imarrayc{1}); colormap(cmap); axis(axstr); colorbar; caxis([mincol maxcol]); daspect(pixdim);

htxt = uicontrol('style','text','String','1','Position',[20 150 30 20]);

%smax = numel(imarrayc);
sstep = [1/(smax-1) 10/(smax-1)];
hslid = uicontrol('style','slider','Min',1,'Max',smax,'Value',1,'Position',...
    [20 200 30 150],'SliderStep',sstep, 'Callback',{@changeimg2, imarrayc, htxt,cmap,mincol, maxcol, pixdim, axstr});
output = 0;
end

function changeimg2(hobj, event,imarray,htxt,cmap,mincol,maxcol, pixdim, axstr)
    ind = int16(get(hobj,'Value'));

    imagesc(imarray{ind});
    axis(axstr); colormap(cmap); colorbar; caxis([mincol maxcol]); daspect(pixdim);

    set(htxt,'String',num2str(ind));
    
end
    