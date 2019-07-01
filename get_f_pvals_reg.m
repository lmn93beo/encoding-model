function [p_vec,F_vec] = get_f_pvals_reg(Xmat, yvec, pred_inds_cell)
% Full model regression
wfull = Xmat\yvec;
error = yvec - Xmat * wfull;
SSE_f = sum(error .^ 2);

n = size(Xmat, 1);
k = size(Xmat, 2);

% Partial regressions
nfactors = length(pred_inds_cell);
for i = 1:nfactors
    % Find partial SSE
    pred_id = pred_inds_cell{i};
    X_part = Xmat;
    X_part(:,pred_id) = [];
    w_part = X_part\yvec;
    error_part = yvec - X_part * w_part;
    SSE_p = sum(error_part .^ 2);
    
    % Calculate F-stat
    m = length(pred_id);
    
    F_vec(i) = (SSE_p - SSE_f) / m / (SSE_f / (n - k));
end

p_vec = nan;