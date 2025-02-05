function cmapLarge = cmap4sense(cmap,cycles,af);


[nc,ny,nx,ns] = size(cmap);

if length(cycles) ~=ns
    disp('Dimensions do not match');
    return
end

nyRed = ns*ny/af;

% shifts need to be integer
shifts = cycles*nyRed/(2*pi);



% intsliceshift = sum(abs(shifts-round(shifts)) > 1e-8) == 0;
%  
% if intsliceshift  
%     shifts = round(shifts);
% else
%     cmapLarge = [];
%     shifts
%     error('shifts are not integer valued, returning')
%     return
% end
% extend FOV by NS
cmapLarge = zeros(nc,ns*ny,nx);

% loop over slices and shift the individual slices according to their 
for k = 1:ns,
    cmapLarge(:,(k-1)*ny+1:k*ny,:) = circshift_(squeeze(cmap(:,:,:,k)),[0 shifts(k) 0]);        
end

% shift to center of large FOV
cmapLarge = circshift(cmapLarge,[0 (ns*ny-ny)/2 0]);
