function A = A_fhp3D(z, nufft_objects)
%
% A = A_fhp3D(z, nufft_objects)
% 
% Defines the forward operator
% 
%   z: the set of images
% 
%   nufft_objects: contains the NUFFT transformations for each image
% 
% Written by Joerg Doepfert 2013
%
[Nimages]=length(nufft_objects);
for i=1:Nimages
         %  do transform for ith image
         A(:,:,i)=nufft_objects{i}*z(:,:,i);
end

end
