function S_inv = build_inverse_covariance(s_mat, nDim)
    % If cov is non-singular, use its inverse.
    % Otherwise, diagonal with 1/cov of each column as in original script.

    if size(s_mat,1) < 2
        % Fallback if too few samples
        S_inv = eye(nDim);
        return;
    end

    C = cov(s_mat);
    if det(C) ~= 0
        S_inv = inv(C);
        return;
    end

    % Diagonal fallback
    S_diag = zeros(nDim,1);
    for idx = 1:min(nDim, size(s_mat,2))
        col = s_mat(:,idx);
        if numel(col) > 1 && var(col) > 0
            S_diag(idx) = 1/var(col);
        else
            S_diag(idx) = 0;
        end
    end
    S_inv = diag(S_diag);
end
