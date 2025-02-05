%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Housekeeping
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ParallelImagingSimulationAndRecon(imgAllSlices,cmapAllSlices,slices,af,caipifactor)
% clear all,
% close all

%%            load images and coilmaps                                   %%

% T1 FLASH example

% imgAllSlices  = importdata('DATA/imgT1fl2d_tra.mat');
% cmapAllSlices = importdata('DATA/cmapT1fl2d_tra.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The matrices have 4 dimensions and following structure
% Nc x Ny x Nx x Ns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% SETUP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Choose slices the acceleration factor and the phace cycles
%
% start for example with ns = 3, af = 3, cycles = [0 0 0]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if caipifactor==0,
    caipifactor=1;
end;

cycles  = mod((0:(length(slices)-1))/caipifactor*2*pi,2*pi);

% Some usefull options
gmax = 6;           % maximum g-factor to display
dispAll = 'on';    % 'on/off' display (not) all figureJs (only recon results),

img  = imgAllSlices (:,:,:,slices);
cmap = cmapAllSlices(:,:,:,slices);

img = img./max(sos(img(:,:)));
cmap = cmap./max(sos(cmap(:,:)));

[nc, ny , nx, ns] = size(img);

if strcmp(dispAll,'on')
    figure(1),set(gcf,'Name', 'Slices and Sensitivity Maps','Color',[1 1 1])
    set(gcf,'pos',[2,378,560,420])
    for k=1:ns,
        %         subplot(2,ns,k),    imagesc(sos(img(:,:,:,k)));  axis off, title(['Slice ' num2str(k)]),axis equal
        %         subplot(2,ns,ns+k), imagesc(sos(cmap(:,:,:,k))); axis off, title(['Sensitivity Slice ' num2str(k)]),axis equal
        
        subplot(2,ns,k),    imab(flipdim(permute(img(:,:,:,k),[3 2 1]),2));  axis off,
        %         title(['Slice ' num2str(k),' could have been acquired with coil array']),axis equal
        title((sprintf('Slice %d  \n acquired \n with coil array',slices(k)))),axis equal
        subplot(2,ns,ns+k), imab(flipdim(permute(cmap(:,:,:,k),[3 2 1]),2)); axis off,
        title(['Sensitivity Slice ' num2str(slices(k))]),axis equal
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The full FoV of one slice is ny
%
% The FoV of the pMSI experiment will be nyRed
% if the acceleration factor equals the number of slices then nyRed = ny
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nyRed  = ns*ny/af;
% adapt parameters if necessary such that nyRed is integer valued
intCheck(nyRed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% These are the relative shifts of the individual slices at the given
% phase cycles in the reduced FOV nyRed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shifts    = cycles*nyRed/(2*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulate the pMSI experiment from fully sampled slice images:
% The script pMSI.m does this job for you. The script applies linear
% phase ramps to the individual slice k-spaces according to the specified
% phase cyles along the specified dimension. In this example along
% the 2nd dimension
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imgMSI = pMSI(img,2,cycles,af);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot the result and repeat with different number of slices,
% acclereation factors and phase cycles.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(dispAll,'on')
    if ns ==1
        figure,
        set(gcf,'Name', 'Simulated accelerated slice','Color',[1 1 1])
    else
        figure,
        set(gcf,'Name', 'Simulated accelerated SMS CAIPIRINHA experiment','Color',[1 1 1])
    end;
    set(gcf,'pos',[520,378,560,420])
    imagesc(sos(imgMSI)); axis off, colormap gray,axis equal
    if ns ==1
        title(sprintf('%d Fourier Transform of the acquired data \n and there is clear wrapping on the phase encoding direction \n where acceleration was performed',ns))
    else
        title(sprintf('%d Slices are overlapped \n and shifted by different amounts \n within the FOV',ns))
    end
end

% PART I: SENSE

% Go back to SETUP and choose data, slices, cycles, accleration and run
% pMSI.m with these parameters to simulate the experiment. Try with and
% without CAIPIRINHA shifts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  pMSI Reconstructions using standard SENSE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% You can use standard inplane SENSE (sense.m) to perform the reconstruction!
% To this end, the sensitivity maps need to be reordered along the phase
% encoding direction within an extended FOV (LargeFOV = Ns*ny) according to
% their shifts in the reduced FOV. The function cmap4sense.m does this job
% for you. Edit the file for details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmapLarge = cmap4sense(cmap,cycles,af);

% imgLarge  = cmap4sense(img,cycles,af);


% If the same proceed to the SENSE recon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Use cmapLarge and af as input for standard SENSE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[reconSENSE, gfactSENSE] = sense(imgMSI,cmapLarge,af);

% Display Recon result and compare with imgLarge

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Take the SENSE results (reconSENSE and gfactor SENSE)
%       crop out the slices and shift them back to the original orientation
%       Hint: Remember the procedure as we generated cmapLarge and imgLarge
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


reconSENSE = circshift(abs(reconSENSE),[-(ns*ny-ny)/2 0]);
gfactSENSE = circshift(abs(gfactSENSE),[-(ns*ny-ny)/2 0]);

% figure(6),
figure,
set(gcf,'Name', ['Parallel imaging (SENSE) recon', ],'Color',[1 1 1])
set(gcf,'pos',[1040,378,560,420])
for k=1:ns,
    subplot(2,ns,k),    imagesc(abs(circshift_(reconSENSE((k-1)*ny+1:k*ny,:),[-shifts(k) 0])),[0 1]);
    axis off,
    if ns==1
        title(sprintf('Accel. factor = %d \n slice = %d' ,(af),slices(k))),
    else
        title(sprintf('Accel. factor = %d \n CAIPI factor = %d \n slice = %d' ,(af),caipifactor,slices(k))),
    end;
    colormap(gray),  axis equal, freezeColors
    subplot(2,ns,ns+k), imagesc(abs(circshift_(gfactSENSE((k-1)*ny+1:k*ny,:),[-shifts(k) 0])),[1 gmax]);
    title('G-factor penalty'), axis off, colormap(jet), axis equal
end
colorbar('East')





