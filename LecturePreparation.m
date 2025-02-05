
%path for the needed files (Matlab and data)

CurrDirectory=pwd;
addpath(genpath(CurrDirectory));



%% Exercise 1 load the single slice k-space data and compute the image 
clear 

% 1a load the data
load kspacedata
close all


% 1b look at the kspace data using imagesc
figureJ(1)
set(gcf,'Position',[0 0 1222 418 ])
subplot(131)
imagesc(abs(kspace)), axis off, axis equal
title('absolute of k-space')

subplot(132)
imagesc(log(abs(kspace)),[2 8]), axis off, axis equal
title('log absolute of k-space')

% 1c compute the image using the matlab function fft and fftshift instead of my function ifft2s

Im=ifft2s(kspace);

figure(1)
subplot(133)
imagesc(abs(Im)), axis off, axis equal
title('absolute of image')
fontScale(1.4)

%% Exercise 2 whoever prescribed the protocol asked for a very large FOV;
% create a reduced field of view by downsampling the k space data 
% specify the indexes of kspace you want to use.
% by how much could you have reduced the the acquisition time if this was a
% 2D or a 3D acquisition?

redFactorx = 4;
redFactory = 5;

kspace_redFOV = kspace(1:redFactorx:end,1:redFactory:end);
% by increasing the spacing between succesive points in k-space by redFactor,
% the FOV is reduced by this same factor

Im_redFOV = ifft2s(kspace_redFOV);
% 
dims    = size(kspace_redFOV);

% alternatively you could simply take the larger FOV image and crop it
% directly in the image space
Im_redFOV_bycropping=crop(Im,dims);
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


% what differences do you see? Why does the data sudenly look so much noisier?
% if you would have changed the FOV simply by croppping the image in the
% real space
colormap gray
%% Exercise 3 given the amount of SNR available in the new FOV dataset, it is maybe better to have a lower resolution currently the prescribed resolution was 0.85... lets make it 1.7
RedRes  = 2;
dimsnew = round(dims/RedRes);

kspace_redFOV_redRES =...
    crop(kspace_redFOV,dimsnew);

% note that by reducing the extent of k-space (done by only taking the
% central part of the k-space), the resolution is automatically reduced.

Im_redFOV_redRES = ifft2s(kspace_redFOV_redRES);


figureJ(3)
set(gcf,'Position',[0 0 1222 418 ])

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
imagesc(abs(Im_redFOV_redRES)), axis off, axis equal, colorbar ('south')
set(gca,'TickLength',[0 0],'YTick',[],'XTick',[],'Color',[1 1 1]) 
axis equal,axis tight, 
ylabel(['reduced by factor ', num2str(redFactorx) ])
xlabel(['reduced by factor ', num2str(redFactory) ])


title(['Image with reduced FOV and red Res by a factor',num2str(RedRes)])

%% Exercise 4 Try to formulate the dependence of SNR in the number of phase encoding  steps used and the resolution;
% Given that the top left square in the image is always background signal, this region will be used to estimate the level of noise of this image.
temp  = (Im( 1:15 , 1:15 ));
Noise = std(temp(:));
temp  = (Im_redFOV( 1:15 , 1:15 ));
Noise_redFOV = std(temp(:));
temp  = (Im_redFOV_redRES( 1:15 , 1:15 ));
Noise_redFOV_redRES = std(temp(:));


display(['the increase of noise level when reducing the number of PE by ', num2str(redFactorx*redFactory)]);
display([' was of ',num2str(Noise_redFOV/Noise)]);

display([' the dependence of the SNR on the number of phase encoding steps is proportional to the square root of the phase encoding steps']);
display([' you can use other reduction factors to see if this relationship holds']);


display(['the decrease of noise level when increasing the volume by ', num2str(RedRes*RedRes)]);
display([' was of ',num2str(Noise_redFOV_redRES/Noise_redFOV)]);
display(['but you should also factor in that you decreased the number of phase encoding steps']);


% as a measure of the signal amplitude you can use the cursor to select a
% region of interest in the tissue of your choice and see how did the
% signal change. (altenatively you can simply look at the range shown by the colorbar)




%%  now that you know everything about k-space lets try to see what is where
% Exercise 5
dimsnew=size(kspace_redFOV_redRES);
mask=zeros(dimsnew);

%Coordinate which information we want to check... you can change the
%coordinate position by changing ks_coord and ky_coord
kx_coord=10;
ky_coord=-10;
% if  n_neighbours = zero, only the info in that coordinate is kept, 
% if a bigger number is the n_neighbours pixels on eaach side are also kept
n_neighbours=0;
mask(kx_coord+round(dimsnew(1)/2+1)+[-n_neighbours:n_neighbours],ky_coord+round(dimsnew(2)/2+1)+[-n_neighbours:n_neighbours])=1;


figureJ(2)
set(gcf,'Position',[ 680   558   340   420])
subplot(221)
imagesc(abs(kspace_redFOV_redRES)), axis off, axis equal
title('absolute of k-space')

subplot(222)
imagesc(abs(kspace_redFOV_redRES.*mask)), axis off, axis equal
title('absolute of k-space')

subplot(223)
imagesc(abs(Im_redFOV_redRES)), axis off, axis equal
title('absolute of image')

% calculates image associated with the masked k-space
Im_masked=ifft2s(kspace_redFOV_redRES.*mask);

subplot(224)
imagesc(real(Im_masked)), axis off, axis equal
title('real part of masked k-space')

colormap(gray)

%%


dimsnew=size(kspace_redFOV_redRES);

%Coordinate which information we want to check... you can change the
%coordinate position by changing ks_coord and ky_coord

k=0
n_neighbours=0;
figureJ(3)
set(gcf,'Position',[ 5    93   666   896]);
for kx_coord=-9:6:9;
   for ky_coord=-9:6:9;
        mask=zeros(dimsnew);
% if  n_neighbours = zero, only the info in that coordinate is kept, 
% if a bigger number is the n_neighbours pixels on eaach side are also kept
mask(kx_coord+round(dimsnew(1)/2+1)+[-n_neighbours:n_neighbours],ky_coord+round(dimsnew(2)/2+1)+[-n_neighbours:n_neighbours])=1;
Im_masked=ifft2s(kspace_redFOV_redRES.*mask);

k=k+1;
set(gcf,'Position',[ 25    70   832   920])

subplot(4,4,k)
imagesc(real(Im_masked)), axis off, axis equal
title(['k = [',num2str(kx_coord),' , ',num2str(ky_coord), '] ' ])

colormap(gray)


   end;
end;


%%
