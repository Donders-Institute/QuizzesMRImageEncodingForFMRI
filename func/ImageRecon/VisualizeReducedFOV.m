function VisualizeReducedFOV(Im,Im_redFOV,Im_redFOV_bycropping,redFactorx,redFactory)

figureJ(1)
set(gcf,'Position',[0 0 1222 418 ],'Name','Exercise 2')

subplot(131)
imagesc(abs(Im)), axis off, axis equal, colorbar('south')
title('Image')

subplot(132)
hold off
imagesc(abs(Im_redFOV)), 
set(gca,'TickLength',[0 0],'YTick',[],'XTick',[],'Color',[1 1 1]) 
axis equal,axis tight, colorbar('south')
title('Image with reduced FOV')
ylabel(['reduced by factor ', num2str(redFactorx) ])
xlabel(['reduced by factor ', num2str(redFactory) ])


subplot(133)
hold off
imagesc(abs(Im_redFOV_bycropping)), 
set(gca,'TickLength',[0 0],'YTick',[],'XTick',[],'Color',[1 1 1]) 
axis equal,axis tight, colorbar('south')
title('Image by cropping edges')
ylabel(['reduced by factor ', num2str(redFactorx) ])
xlabel(['reduced by factor ', num2str(redFactory) ])
colormap gray
