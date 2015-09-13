function X_k = SVT_LS(Y, A, At, tau, delta, lambda, k_max, tol, mask)
%
% X = SVT_LS(Y, A, At, tau, delta, lambda, k_max, tol, mask)
% 
% SVT-LS reconstruction algorithm as outlined in the main article
% 
% 
% 
% Code written by Joerg Doepfert 2013, based on SVT.m by Emmanuel Candes 
% and Stephen Becker, http://svt.stanford.edu/code.html


    % initialization    
    Y_kminus1 = zeros(size(Y));           % Y_0
    M = repmat(mask, [1, 1, size(Y,3)]);  % signal mask
    normY = norm(Y(:));

   for k = 1:k_max,

            X_kminus1 = At(Y_kminus1);

            % LS constraint
            Z_k = X_kminus1 .* M;
            X_k = X_kminus1 + lambda * (Z_k-X_kminus1);

            % singular value shrinkage
            [X_k, r, sigma] = shrink(X_k, tau);


            % data fidelity
            A_X_k = A(X_k);
            Y_k = Y_kminus1 + delta * (Y-A_X_k);

            % update Y_kminus1 for next iteration
            Y_kminus1 = Y_k;

            % stop if tolerance level is reached
            relRes = norm(A_X_k(:) - Y(:)) / normY;
            if (relRes < tol)
                break
            end

            % stop if algorithm diverges
            if (norm(A_X_k(:) - Y(:)) / normY > 1e5)
                disp('Divergence!');
                break
            end

            % plot reco progress
            residual(k) = relRes;
            rank(k) = r;
            % fprintf('iter %4d, rank %2d, rel. residual %.3e\n', k, r, relRes);

            subplot(1,2,1);
            ShowImages(abs(X_k));
            title(['current guess X_k, k = ' num2str(k)]);
            
            subplot(1,2,2);
            [AX,H1,H2] = plotyy(1:k, residual(1:k), 1:k, rank(1:k));
            set(get(AX(1),'Ylabel'),'String','relative residual') 
            set(get(AX(2),'Ylabel'),'String','rank of X_k') 
            xlabel('k')
            title(['rel residuum / rank: ',num2str(residual(k)),' / ',num2str(rank(k))]);
            drawnow;
    end

end


