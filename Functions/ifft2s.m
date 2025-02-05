function out = ifft2s(data)
% function out = ifft2s(data)
% 
% Do 2D FFT on each slice of data and fftshift

out = fftshift(fftshift(ifft2(ifftshift(ifftshift(data,1),2)),1),2);
