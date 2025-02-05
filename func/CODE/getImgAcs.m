function imgAcs = getImgAcs (img,nacs);


[nc, ny, nx, ns ] = size(img);

if nargin < 2
    nacs(1) = 32;
    nacs(2) = 32;
end

disp (['ACS size will be ' num2str(nacs(1)) ' x ' num2str(nacs(2)) ])

 sig = fftshift(fftshift(fft(fft(fftshift(fftshift(img,2),3),[],2),[],3),2),3);
idxy = (ny-nacs(1))/2 + [1:nacs(1)];
idxx = (nx-nacs(2))/2 + [1:nacs(2)];
sigAcs = sig(:,idxy,idxx,:);
imgAcs = ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(sigAcs,2),3),[],2),[],3),2),3);

