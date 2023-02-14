function [ERROR_CALC,RESULTS] = getErrorCalculations(SIM,FLAG,N,P,RESULTS,ESTIMATOR)
% Initialize Variables
    plant_states = RESULTS.Slink_plant.xNR.x_soln;
    est_states   = RESULTS.EST.ASY.x_soln;
    
    [N_states, N_samples] = size(plant_states);
    
    plant_outputs = RESULTS.Slink_plant.zNA.z_soln;
    est_outputs   = RESULTS.EST.ASY.z_soln_ALL;
    
    [N_outputs, ~] = size(plant_outputs);
    
    
    idx = floor(N_samples/2);

% Calculate the suspected error in the states
    error_states = est_states(:,idx:end) - plant_states(:,idx:end);
    
    covar_states = nan(N_states);
    for i = 1:N_states
        for j = 1:N_states
            [covar_states(i,j)] = calcPopCovar(error_states(i,:), error_states(j,:));
        end
    end
    
% Calculate the error in the outputs
    error_outputs = est_outputs(:,idx:end) - plant_outputs(:,idx:end);
    
    covar_outputs = nan(N_outputs);
    for i = 1:N_outputs
        for j = 1:N_outputs
            [covar_outputs(i,j)] = calcPopCovar(error_outputs(i,:), error_outputs(j,:));
        end
    end



% Compare Calculated to Expected
    P_calc = covar_states;
    P_exp = ESTIMATOR.P_infty;

%     disp(newline)
%     disp('Cov(states) calc is')
%     disp(num2str(covar_states))
%     disp('P_\infty is ')
%     disp(num2str(P_exp))
%     disp('With difference of')
%     disp(num2str(P_calc-P_exp))
%     disp('Ratio of error')
%     disp(num2str(abs((P_calc-P_exp)./P_exp)))

    CPCT_calc = covar_outputs;
    CPCT_exp = ESTIMATOR.C_des_ROM * ESTIMATOR.P_infty * ESTIMATOR.C_des_ROM';

%     disp(newline)
%     disp('Cov(outputs) calc is')
%     disp(num2str(covar_outputs))
%     disp("CP_\inftyC' is ")
%     disp(num2str(CPCT_exp))
%     disp('With difference of')
%     disp(num2str(CPCT_calc-CPCT_exp))
%     disp('Ratio of error')
%     disp(num2str(abs((CPCT_calc-CPCT_exp)./CPCT_exp)))


    for i = 1:3
        CPCT_exp_diag(i) = ESTIMATOR.C_des_ROM(i,:) * ESTIMATOR.P_infty * ESTIMATOR.C_des_ROM(i,:)';
    end
%     disp(newline)
%     disp(newline)
%     disp(num2str(CPCT_exp_diag))

    ERROR_CALC.P_calc = P_calc;
    ERROR_CALC.P_infty = P_exp;
    ERROR_CALC.CPC_calc = CPCT_calc;
    ERROR_CALC.CPC_infty = CPCT_exp;
    ERROR_CALC.CPC_infty_diag = CPCT_exp_diag;


end

function [mu] = calcMean(x)
    mu = mean(x);
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
