function out = sos(data,dim);

% Calculates "Sum-Of-Squares"-Image (SOS) along dimension dim
% Call: [SOS] = SOS(DATA,dim)
% default dim is 1 
%
% Input:        Multidemensional data set 'DATA' with coil dimension at dim
% Output:      "Sum-Of-Squares"-Image 'SOS'

if nargin <2, 
    dim=1; 
end

out=squeeze(sqrt(mean(abs(data).^2,dim)));
