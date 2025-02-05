function imgMSI = pMSI(img_in,dimCycle,cycles,af)


% This code applies linear physe cycles to individual slices and allows for the simulation of Multiband CAIPIRINHA acquisitions  
%
% INPUT:  img_in                        can be of any size: Slice dimension should be last
%                                       eg. [Ncoils x Ny x Nx x Ns]
%         dim                           specifies the dimeinsion on which
%                                       the phase cycle is applied
%         cycle                         cycle should be a vector of length Ns 
%                                       eg. [dphi1 dPhi2 ...]
%
%         af                            integer acceleration factor
%                                       eg. af = 1,2,3,4 ...
%                                       Ns*Ny must be divideable by af
%
%
% Example 1:  img_CAIPI = pMSI(img,2,[0 pi/2],2)
%
%
%
%
%
% Example 2:
%
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check the Input data for correct settings
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ns = length(cycles);
Nd = ndims(img_in);

if (size(img_in,Nd) ~= Ns && Ns > 1 || Nd < 2)
    disp(['Check Input...Number of slices and sycles are not the same']);
    imgMSI = [];
    return 
end

dimSlice = Nd;

if (Ns*size(img_in,dimCycle)/af ~= round(Ns*size(img_in,dimCycle)/af))
    disp(['Check acceleration factor ... Ny*Ns/af should be integer valued']);
    imgMSI= [];
    return 
end

if (Ns == 1)
    dimSlice      = dimSlice +1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reorder matrix: Cycle are applied to 1st Slices are 2nd
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




[img_in, order] = reorder(img_in,dimCycle,dimSlice);

sz = size(img_in);
Ny = sz(1);
Ns = sz(2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create LargeFOV image and put slices in center
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imgLargeFOV = repmat(zeros(size(img_in(:,:,:))),[Ns 1 1]);
imgLargeFOV(Ns*Ny/2-Ny/2+1:Ns*Ny/2+Ny/2,:,:) = img_in(:,:,:);

clear img_in;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fourier transformation along first dimension
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigLargeFOV = fftshift(fft(fftshift(imgLargeFOV,1),[],1),1);

clear imgLargeFOV;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Apply acceleration  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sigRedFOV = sigLargeFOV(1:af:end,:,:);

clear sigLargeFOV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Apply Phasecycles to the individual slices: can be done much faster
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigMSI = zeros(Ny*Ns/af,prod(sz(3:end)));

% for m = 1:Ny*Ns/af;
%     
%     for k = 1:Ns,
%          sigTmp(k,1,:) = sigRedFOV(m,k,:).*exp(-1j*(m-1)*cycles(k));
%     end
%     
%     sigMSI(m,:) =  sum(sigTmp,1); % Take complex sum along slices
% 
% end
% sigTmp = bsxfun(@times,bsxfun(@times,permute((1:Ny*Ns/af),[2,1]),1:Ns),sigRedFOV);
 sigMSI =  squeeze(sum(bsxfun(@times,exp(-1j*bsxfun(@times,permute((1:Ny*Ns/af),[2,1]),cycles(1:Ns))),sigRedFOV),2)); % Take complex sum along slices


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inverse Fourier transformation along first dimension
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imgMSI = ifftshift(ifft(ifftshift(sigMSI,1),[],1),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reshape and reorder 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imgMSI = reshape(imgMSI,[Ns*Ny/af 1 sz(3:end)]);
imgMSI = permute(imgMSI,order);



imgMSI(abs(imgMSI)<0.001*max(abs(imgMSI(:))))=0;



end

%%




%%

function [img_out,order] = reorder(img_in,dimCycle,dimSlice);

Nd      = ndims(img_in);        % Number of dimensions of input data

if dimSlice > Nd;
    Nd = dimSlice;
end

perm(1)   = dimCycle;
perm(2)   = dimSlice;

if dimCycle > 1
  perm(2+1:2+dimCycle-1)   = 1:dimCycle-1;
end

perm(2+dimCycle:Nd)       = dimCycle+1:Nd-1;

img_out = permute(img_in,perm);

[~ , order] = sort(perm);

end


