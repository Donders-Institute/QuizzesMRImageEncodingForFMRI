function [recon,gfact] = sense(imfold,cmap,af)

%   SENSE recon for 1d acceleration
%
%   IN:         img              folded coil images         (#coils, ny/af, nx)
%               cmap             coilmaps                   (#coils, ny, nx)    
%                 af             Acceleration factor        (integer)
%               
%
%   OUT:       recon              Reconstructed images      (ny,nx)    
%              gfact              Geometry factor           (ny,nx)
%
%   If you use this code please reference Pruessmann et al;



%warning('on','last');


cmap   = circshift(cmap,[0 -size(cmap,2)/2  0]);
imfold = circshift(imfold  ,[0 -size(imfold,2)/2 0]); 

[nc,ny,nx]    = size(cmap);
[ncr,nyr,nxr] = size(imfold);

if nx~=nxr || nyr*af~=ny || nc~=ncr
    disp('Dimensions of imfold and cmap do not match!')
    return;
end

cmap(isinf(cmap))=0;        %Make sure there are no bad entries in the coil maps
cmap(isnan(cmap))=0;


recon = zeros(ny,nx);
gfact = zeros(ny,nx);


lambda = 1e-3;


for x=1:nx
    
    if round(x/(nx/4))== x/(nx/4), fprintf('...%d', x),end
    
    for y=1:nyr
        % sensitivity matrix
        c = squeeze(cmap(:,y:nyr:ny,x));
        sc = max(abs(c(:)));
        % invert
         cc = c'*c;
        
        % regularize
        s = svd(cc,0);
        s = sqrt(max(abs(s)));
        
        if s > 0.1
        u = inv(cc+eye(size(cc)).*lambda.*s.^2); 
        
        % unfolding matrix
        unfold = u*c';   %
        
%         unfold = pinv(c);
                
        % unfold aliased pixels      
        r= unfold*imfold(:,y,x);
        % calculate g-factor (Pruessmann et al)
        g=sqrt(diag(cc).*diag(u)); 
   
        % populate image and g-factor
        recon (y:nyr:ny,x)=r;
        gfact (y:nyr:ny,x)=g;
         end
    end
end


gfact = circshift(gfact,[size(gfact,1)/2 0]);
recon = circshift(recon,[size(recon,1)/2 0]);

fprintf('\n');

end



