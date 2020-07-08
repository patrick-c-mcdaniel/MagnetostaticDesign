function [] = displayStreamFunctionAadj(triangle,s,node)

figure

%First approximation based on the mean on each triangle
TriAmpl = zeros(size(triangle,1),1);
for i=1:size(triangle,1)
    TriAmpl(i) = (s(triangle(i,1))+s(triangle(i,2))+s(triangle(i,3)))/3;
end
AREAs = zeros(size( triangle,1 ),1);
for iarea = 1:size( triangle,1 )
    AREAs(iarea) = 1/2*norm( cross( node(triangle(iarea,2),:)-node(triangle(iarea,1),:), ...
                                    node(triangle(iarea,3),:)-node(triangle(iarea,1),:) ),2 );
    
end


trisurf(triangle,node(:,1),node(:,2),node(:,3),TriAmpl./((AREAs).^(1/3)))
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');
axis equal;
%caxis auto
%caxis([-2e5 2e5])

v = caxis;
caxis([-max(abs(v)) max(abs(v))]);
cmap = colormap;
newColorMap = gaelColormap(cmap);
colormap(newColorMap);