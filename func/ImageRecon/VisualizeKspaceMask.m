function VisualizeKspaceMask(kspace_redFOV_redRES , mask ,Im_redFOV_redRES  , Im_masked)

% close(2)
figureJ()
subplot(221)
imagesc(abs(kspace_redFOV_redRES)), axis off, axis square
title('absolute of k-space')

subplot(222)
imagesc(abs(kspace_redFOV_redRES.*mask)), axis off, axis square
% imagesc(abs(kspace_redFOV_redRES.*mask)), axis square
title('absolute of k-space')
xlabel('ky')
ylabel('kx')

subplot(223)
imagesc(abs(Im_redFOV_redRES)), axis off, axis equal
title('absolute of image')


subplot(224)
imagesc(real(Im_masked)), axis off, axis equal
xlabel('y')
ylabel('x')
title('real part of masked k-space')

colormap(gray)