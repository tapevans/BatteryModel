%% Calculated Covariance

function [Pcalc] = getP_calc(y_plant,y_est,idx)
% Covariance
        idx = 1000;
        [N_states, ~] = size(y_plant);

        
        y_plant = y_plant';
        y_est   = y_est';
        % Error Calc (Deviation)
        for i = 1:N_states
            error(:,i) = y_est(idx:end,i) - y_plant(idx:end,i);
        end
        error = error';
        
        % Covar calc
        covar = nan(N_states);
        for i = 1:N_states
            for j = 1:N_states
                [covar(i,j)] = calcPopCovar(error(i,:), error(j,:));
            end
        end

        Pcalc = covar;

end


function [Covar] = calcPopCovar(error_x, error_y)
    Covar = (error_x*error_y')/length(error_x);
end