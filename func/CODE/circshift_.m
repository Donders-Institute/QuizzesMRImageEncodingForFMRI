
function Y = circshift_(X,V,dim)

%   Same as circshift but allows for non-integer shifts
%   Y = CIRCSHIFT(X,V) circularly shifts the values in the array X
%   by V elements. V is a vector of integers where the N-th element
%   specifies the shift amount along the N-th dimension of
%   array X.



if nargin < 3

    if length(V) ~= ndims(X);
        error('Check Input');
        return 
    end
    
    
    for k = 1:ndims(X);
        X = nonIntCircshift_(X,V(k),k);
        
    end
    
    Y=X;
    
else
    
    if length(V) > 1 || dim > length(V) || dim < 1;
        error('Check Input');
        return 
    end
    
    Y = nonIntCircshift_(X,V,dim);
    
end

end

function Y = nonIntCircshift_(X,V,dim)
sz = size(X);
sz(dim) = 1;
rs = ones(1,ndims(X));
rs(dim) = size(X,dim);
sig = fftshift(fft(fftshift(X,dim),[],dim),dim);
phase = repmat(reshape(exp(-1i*2*pi/size(X,dim)*V.*[1:size(X,dim)]),rs),sz);
sig = sig.*phase;
Y = ifftshift(ifft(ifftshift(sig,dim),[],dim),dim);
end