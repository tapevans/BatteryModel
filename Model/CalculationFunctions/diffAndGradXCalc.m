%% diffAndGradCalc
function [X_diff, X_grad] = diffAndGradXCalc(X, SIM__diff_CV_x_vec_inv)
    diff   = X(2:end) - X(1:end-1);
    X_diff = [nan diff nan];
    X_grad = X_diff.*SIM__diff_CV_x_vec_inv;
end