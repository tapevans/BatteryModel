%% Data2Covariance
function [Covar] = data2Covar(data)
% Initialize Outputs
    [ N_outputs , ~] = size(data); % Number of Outputs
    Covar = nan(N_outputs);

% Calculate the Mean 
    [mu] = calcMean(data); %??? will mu be a vector

% Calculate the Error
    [error] = calcError(data,mu); %??? will this handle a vector

% Calculate Covariance
    for i = 1:N_outputs
        for j = 1:N_outputs
            [Covar(i,j)] = calcPopCovar(error(i,:), error(j,:));
        end
    end
end


%% Helper Functions
function [mu] = calcMean(x)
    mu = mean(x,2); % ,2 is to take the mean across the rows and not the columns
end

function [error] = calcError(x,mu)
% difference between x and the expected value (mu)
    error = x-mu;
end

function [Covar] = calcPopCovar(error_x, error_y)
    Covar = (error_x*error_y')/length(error_x);
end

function [Covar] = calcSampleCovar(error_x, error_y)
    Covar = (error_x*error_y')/(length(error_x)-1);
end