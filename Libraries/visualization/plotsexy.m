function [phn] = plotsexy(varargin)
    
    phn = plot(varargin{:});
   % htit = title('Magnitude');
   % set(gca, 'xdir', 'reverse')
    set(phn , 'MarkerSize' , 10,'LineWidth',2);
    %set(hj2, 'Color',[0 0 0.7],'LineWidth',2);
   % set(htit, 'FontName','Liberation Sans','FontSize',18); 
    set(gca,'FontName','Liberation Sans','FontSize',14,'LineWidth',2);
    set(gcf,'Color','white')
    grid on; grid minor
   % ylim([0 5]);
end