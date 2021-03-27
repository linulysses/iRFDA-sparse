function M = lc_mean(V1,V2,V3,L1,L2,L3,mask)

    dim = size(mask);
    N = sum(sum(sum(mask)));
    S = zeros(3,3,N);
    k = 0;
    for x = 1:dim(1)
        for y = 1:dim(2)
            for z = 1:dim(3)
                if mask(x,y,z) == 1 %&& min([L1(x,y,z),L2(x,y,z),L3(x,y,z)]) > 1e-6
                    v1 = squeeze(V1(x,y,z,:));
                    v2 = squeeze(V2(x,y,z,:));
                    v3 = squeeze(V3(x,y,z,:));
                    k = k + 1;
                    S(:,:,k) =  L1(x,y,z) * (v1 * v1') + L2(x,y,z) * (v2 * v2') + L3(x,y,z) * (v3 * v3');
                end
            end
        end
    end    
    %S = S(:,:,1:k);
    M = frechet_mean_LogCholesky(S);
    
%% helper functions

% The Cholesky factor (lower triangular) of an SPD
% @param S a SPD
    function L = ch(S)
        L = chol(S)';
    end

    function M = frechet_mean_LogCholesky(S)
        R = 0;
        D = 0;
    
        n = size(S,3);
        m = 0;
        for i = 1:n

            if ~(any(any(isinf(S(:,:,i)))) || min(eig(S(:,:,i))) <= 0)
                C = chol(S(:,:,i))'; %ch(S(:,:,i));
                R = R + tril(C,-1);
                D = D + log(diag(C));
                m = m + 1;
            end
        end
        R = R / m;
        D = D / m;
        L = R + diag(exp(D));
        
        M = L * L';
    
    end

end