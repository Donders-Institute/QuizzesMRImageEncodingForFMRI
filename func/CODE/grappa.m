function [recon , gfact ] = grappa(imfold,imacs,af);


imsrc = imacs;
imtrg = imacs;


[nc    nyr   nx   ]   = size(imfold);
[ncsrc nysrc nxsrc]   = size(imsrc);
[nctrg nytrg nxtrg]   = size(imtrg);

if nc~=ncsrc || nysrc~=nytrg || nysrc~=nytrg
    disp('Dimensions of imfold and cmap do not match!')
    return;
end

ny = nyr*af;

% Fourier transform in k-space
src = fftshift(fftshift(fft(fft(fftshift(fftshift(imsrc,2),3),[],2),[],3),2),3);
trg = fftshift(fftshift(fft(fft(fftshift(fftshift(imtrg,2),3),[],2),[],3),2),3);


kernel  = calcKernel    (src,trg,af);
weights = getWeights    (kernel,ny,nx);
recon   = applyWeights  (imfold,weights);
gfact   = calcGfact     (weights,imtrg,af);

recon = squeeze(sqrt(sum(abs(recon).^2,1)));

end




%%  Helper functions



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Calculation of GRAPPA weights in kspace from ACS data
%  see Griswold et al
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ws_kernel = calcKernel(acs_src,acs_trg,af)

[nc_src,ny_src,nx_src]=size(acs_src);

nc    = nc_src;
nyacs = ny_src;
nxacs = nx_src;

srcx = 5;                         % should be odd
srcy = 4;                         % can now be even and odd

src=zeros(nc*srcy*srcx,(nyacs-(srcy-1)*af)*(nxacs-(srcx-1)));
trg=zeros(nc*af,(nyacs-(srcy-1)*af)*(nxacs-(srcx-1)));


cnt = 0;  % This is a lazy counter. could be done much faster.

for xind=floor(srcx./2)+1:nxacs-floor(srcx./2),
    for yind=1:nyacs-(srcy-1)*af,
        cnt=cnt+1;
        src(:,cnt)=reshape(acs_src(:,yind:af:yind+(srcy-1)*af,xind-floor(srcx./2):xind+floor(srcx./2)),nc*srcy*srcx,1);
        trg(:,cnt)=reshape(acs_trg(:,yind+floor((af*(srcy-1)+1)/2) - floor(af/2) : yind + floor((af*(srcy-1)+1)/2) - floor(af/2) + af -1 ,xind),nc*af,1);
    end
end


ws = trg*pinv_reg(src,1e-4);


ws_tmp = reshape(ws,[nc,af,nc,srcy,srcx]);                                 %Reshape weight set
% ws_tmp = flip(flip(ws_tmp,4),5);                                     %flip source points in ky and kx for the convolution
ws_tmp = flipdim(flipdim(ws_tmp,4),5);                                     %flip source points in ky and kx for the convolution


for k=1:af,
    ws_kernel(:,:,k:af:af*srcy,:) = ws_tmp(:,k,:,:,:);                 %reconstruction kernel
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                   
%                                                                    
%  Calculate grappa weights for reconstruction in image space                                                                  
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ws_img = getWeights(ws_kernel,ny,nx)

[nc,~,nky,nkx] = size(ws_kernel);

ws_k = zeros(nc,nc,ny,nx);
ws_k(:,:,ceil((ny-nky)/2)+1:ceil((ny+nky)/2),ceil((nx-nkx)/2+1):ceil((nx+nkx)/2)) = ws_kernel;  %put reconstruction kernel in the center of matrix

tmp0 = ifftshift(ws_k,3);           % shift in phase
tmp1 = ifftshift(tmp0,4);           % shift in read
tmp0 = ifft(tmp1,[],3);             % ifft in phase
tmp1 = ifft(tmp0,[],4);             % ifft in read
tmp0 = ifftshift(tmp1,3);           % shift in phase
tmp1 = ifftshift(tmp0,4);           % shift in read
ws_img = ny*nx*tmp1;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Weights application in image space
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function recon = applyWeights(img_red,ws_img);


nc = size(ws_img,1);
ny = size(ws_img,3);
nx = size(ws_img,4);

nyr = size(img_red,2);

if (ny > nyr);
  
    af = round(ny/nyr);
    
    disp('Assuming the data is passed without zeros at not sampled lines ......')
    disp(['Acceleration factor is af = ' num2str(af) '.....'] );
       
    sig_red = fftshift(fftshift(fft(fft(fftshift(fftshift(img_red,2),3),[],2),[],3),2),3);

    sig_new = zeros(nc,ny,nx);
    sig_new(:,1:af:end,:) = sig_red;
    sig_red = sig_new;
    
    img_red = ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(sig_red,2),3),[],2),[],3),2),3);
        
end

recon = zeros(nc,ny,nx);

for k = 1:nc,
    recon(k,:,:) = sum(squeeze(ws_img(k,:,:,:)).*img_red,1);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gfactor calculation according to Breuer et al
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function g = calcGfact(ws_img,imgref,af,R)


nc = size(ws_img,1);
ny = size(ws_img,3);
nx = size(ws_img,4);

ws_img = ws_img/af;

if nargin < 4;
    R = eye(nc);
    disp('No noise correlations has been passed ......')
    disp('Assuming perfect coils without correlations......')
end

disp('Calculating G-factor......')

sigref = fftshift(fftshift(fft(fft(fftshift(fftshift(imgref,2),3),[],2),[],3),2),3);

[nc, nyref, nxref] = size(sigref);

% filter kspace
sigref = bsxfun(@times,sigref,reshape(tukeywin(nyref,1),[1 nyref]));
sigref = bsxfun(@times,sigref,reshape(tukeywin(nxref,1),[1 1 nxref]));
% 
sigreff = zeros(nc,ny,nx);
yidx = floor((ny-nyref)/2) + [1:nyref];
xidx = floor((nx-nxref)/2) + [1:nxref];
sigreff(:,yidx,xidx) = sigref;

imgref = ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(sigreff,2),3),[],2),[],3),2),3);

g = zeros(ny,nx);


Z = eye(nc);

for y = 1:ny,
    for x = 1:nx,
        W = ws_img(:,:,y,x);            % Weights in image space
        tmp = imgref(:,y,x);
        n = tmp'./sqrt(sum(abs(tmp).^2,1));
        g(y,x) = sqrt(abs((n*W)*R*(n*W)'))./sqrt(abs((n*Z)*R*(n*Z)'));       % This is the generalized g-factor formulation
    end
    
end


end


function X = pinv_reg(A,lambda)
%PINV   Pseudoinverse.
%   X = PINV(A) produces a matrix X of the same dimensions
%   as A' so that A*X*A = A, X*A*X = X and A*X and X*A
%   are Hermitian. The computation is based on SVD(A) and any
%   singular values less than a tolerance are treated as zero.
%   The default tolerance is MAX(SIZE(A)) * NORM(A) * EPS.
%
%   PINV(A,TOL) uses the tolerance TOL instead of the default.
%
%   See also RANK.

%   Copyright 1984-2001 The MathWorks, Inc. 
%   $Revision: 5.11 $  $Date: 2001/04/15 12:01:37 $
[m,n] = size(A);

if n > m
   X = pinv_reg(A',lambda)';
else
    
   AA = A'*A;

   S = svd(AA,0);
   
   S = sqrt(max(abs(S)));
   X = (AA+eye(size(AA)).*lambda.*S.^2)\A'; 
  
end

end
