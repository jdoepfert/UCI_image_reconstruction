function At = At_fhp3D(z, nufft_objects) 
%
% At = At_fhp3D(z, nufft_objects, Npix, Nspokes)
% 
% Defines the adjoint of the forward operator
% 
%   z: k-space data
% 
%   nufft_objects: contains the NUFFT transformations for each image
% 
% Written by Joerg Doepfert 2013
%
[Nimages]=length(nufft_objects);
for i=1:Nimages
    %do transform for ith image
    At(:,:,i)=nufft_objects{i}'*z(:,:,i);
end

end

