function [X_new, r, sigma] = shrink( X, tau )
%
% X_new = shrink( X, tau )
% 
% Performs singular value shrinkage on image set X
% 
%   tau  : parameter controlling the amount of shrinkage
%   r    : rank
%   sigma: singular values
% 
% Code written by Joerg Doepfert 2013, based on SVT.m by Emmanuel Candes 
% and Stephen Becker, http://svt.stanford.edu/code.html

   [N1, N2, Nimages] = size(X);
   X_temp = reshape(X, N1 * N2, Nimages);  
   [U, Sigma, V] = svd(X_temp, 'econ');
   sigma = diag(Sigma); r = sum(sigma > tau); 
   U = U(:,1:r); V = V(:,1:r); sigma = sigma(1:r) - tau; 
   Sigma = diag(sigma);   
   X_new = U * Sigma * V';
   X_new = reshape(X_new, N1, N2, Nimages);

end


