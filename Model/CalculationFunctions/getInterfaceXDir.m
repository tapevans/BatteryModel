%% getInterfaceXDir
function [X_interface] = getInterfaceXDir(X, SIM__interp_x_interface)
    diff        = X(2:end) - X(1:end-1);
    shifted     = [nan diff nan].*SIM__interp_x_interface;
    X_interface = shifted + [nan X(1:end-1) nan];
end