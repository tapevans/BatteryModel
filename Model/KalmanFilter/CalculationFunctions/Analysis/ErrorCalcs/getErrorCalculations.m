function [ERROR_CALC,RESULTS] = getErrorCalculations(SIM,FLAG,N,P,RESULTS,ESTIMATOR)
idx = floor((1/2)*FLAG.N_samples)+2;
idx = 500; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Fix this !!!!!!!!!!!!!!!!!!!!!!!!!!!
% idx = floor((3/4)*FLAG.N_samples);

% Initialize Variables
%     plant_states = RESULTS.Slink_plant.xNR.x_soln;
%     est_states   = RESULTS.EST.ASY.x_soln;
%     
%     [N_states, N_samples] = size(plant_states);
      
%     idx = floor(N_samples/2);

% % Calculate the suspected error in the states
%     error_states = est_states(:,idx:end) - plant_states(:,idx:end);
%     
%     covar_states = nan(N_states);
%     for i = 1:N_states
%         for j = 1:N_states
%             [covar_states(i,j)] = calcPopCovar(error_states(i,:), error_states(j,:));
%         end
%     end

% % Compare Calculated to Expected
%     P_calc = covar_states;
%     P_exp = ESTIMATOR.P_infty;
% 
%     if FLAG.Analysis.dispResults
%         disp(newline)
%         disp('Cov(states) calc is')
%         disp(num2str(covar_states))
%         disp('P_\infty is ')
%         disp(num2str(P_exp))
%         disp('With difference of')
%         disp(num2str(P_calc-P_exp))
%         disp('Ratio of error')
%         disp(num2str(abs((P_calc-P_exp)./P_exp)))
%     end
    
% ERROR_CALC.P_calc = P_calc;
% ERROR_CALC.P_infty = P_exp;

%% Calculate the error in the outputs (ASY)
    plant_outputs = RESULTS.EST.PLANT.z_soln_ALL;
    est_outputs   = RESULTS.EST.ASY.z_soln_ALL;
    
    [N_outputs, ~] = size(plant_outputs);
    %error_outputs = est_outputs(:,idx:end)   - plant_outputs(:,idx:end);
    error_outputs = est_outputs(:,idx+1:end) - plant_outputs(:,idx:end-1);
    %error_outputs = est_outputs(:,idx+2:end) - plant_outputs(:,idx:end-2);
    
    covar_outputs = nan(N_outputs);
    for i = 1:N_outputs
        for j = 1:N_outputs
            [covar_outputs(i,j)] = calcPopCovar(error_outputs(i,:), error_outputs(j,:));
        end
    end


    CPCT_calc = covar_outputs;
    if FLAG.EstimatorModel == 1 || FLAG.EST.SepHK == 0
        CPCT_exp = ESTIMATOR.C_des_ROM{1} * ESTIMATOR.P_infty{1} * ESTIMATOR.C_des_ROM{1}';
    else
        for OO = 1:length(ESTIMATOR.P_infty)
            CPCT_exp_temp = ESTIMATOR.C_des_ROM{OO} * ESTIMATOR.P_infty{OO} * ESTIMATOR.C_des_ROM{OO}';
            CPCT_exp(OO,OO) = CPCT_exp_temp(2,2);
        end
    end


    if FLAG.Analysis.dispResults
%         disp(newline)
%         disp('Cov(outputs) calc is')
%         disp(num2str(covar_outputs))
%         disp("CP_\inftyC' is ")
%         disp(num2str(CPCT_exp))
%         disp('With difference of')
%         disp(num2str(CPCT_calc-CPCT_exp))
%         disp('Ratio of error')
%         disp(num2str(abs((CPCT_calc-CPCT_exp)./CPCT_exp)))
    end


%     for i = 1:N_outputs
%         CPCT_exp_diag(i) = ESTIMATOR.C_des_ROM(i,:) * ESTIMATOR.P_infty * ESTIMATOR.C_des_ROM(i,:)';
%     end
    CPCT_exp_diag  = diag(CPCT_exp)';
    CPCT_calc_diag = diag(CPCT_calc)';

    if FLAG.Analysis.dispResults
        disp(newline)
        disp(newline)
        disp('CPC^T diagonal Estimator')
        disp(num2str(CPCT_exp_diag))
        disp(newline)
        disp('CPC^T diagonal Calculated')
        disp(num2str(CPCT_calc_diag))
    end

    
    ERROR_CALC.CPC_calc_asy = CPCT_calc;
    ERROR_CALC.CPC_calc_diag_asy = CPCT_calc_diag;
    ERROR_CALC.CPC_infty_asy = CPCT_exp;
    ERROR_CALC.CPC_infty_diag_asy = CPCT_exp_diag;


%% Calculate the error in the outputs (VAR)
    plant_outputs = RESULTS.EST.PLANT.z_soln_ALL;
    est_outputs   = RESULTS.EST.VAR.z_soln_ALL;
    
    [N_outputs, ~] = size(plant_outputs);
    %error_outputs = est_outputs(:,idx:end)   - plant_outputs(:,idx:end);
    error_outputs = est_outputs(:,idx+1:end) - plant_outputs(:,idx:end-1);
    %error_outputs = est_outputs(:,idx+2:end) - plant_outputs(:,idx:end-2);
    
    covar_outputs = nan(N_outputs);
    for i = 1:N_outputs
        for j = 1:N_outputs
            [covar_outputs(i,j)] = calcPopCovar(error_outputs(i,:), error_outputs(j,:));
        end
    end

    CPCT_calc = covar_outputs;
    if FLAG.EstimatorModel == 1 || FLAG.EST.SepHK == 0
        CPCT_exp = ESTIMATOR.C_des_ROM{1} * ESTIMATOR.P_infty{1} * ESTIMATOR.C_des_ROM{1}';
    else
        for OO = 1:length(ESTIMATOR.P_infty)
            CPCT_exp_temp = ESTIMATOR.C_des_ROM{OO} * ESTIMATOR.P_infty{OO} * ESTIMATOR.C_des_ROM{OO}';
            CPCT_exp(OO,OO) = CPCT_exp_temp(2,2);
        end
    end

    if FLAG.Analysis.dispResults
        disp(newline)
        disp(newline)
        disp('Variable Estimator')
%         disp('Cov(outputs) calc is')
%         disp(num2str(covar_outputs))
%         disp("CP_\inftyC' is ")
%         disp(num2str(CPCT_exp))
%         disp('With difference of')
%         disp(num2str(CPCT_calc-CPCT_exp))
%         disp('Ratio of error')
%         disp(num2str(abs((CPCT_calc-CPCT_exp)./CPCT_exp)))
    end


%     for i = 1:N_outputs
%         CPCT_exp_diag(i) = ESTIMATOR.C_des_ROM(i,:) * ESTIMATOR.P_infty * ESTIMATOR.C_des_ROM(i,:)';
%     end
    CPCT_exp_diag  = diag(CPCT_exp)';
    CPCT_calc_diag = diag(CPCT_calc)';

    if FLAG.Analysis.dispResults
        disp(newline)
        disp(newline)
        disp('CPC^T diagonal Estimator')
        disp(num2str(CPCT_exp_diag))
        disp(newline)
        disp('CPC^T diagonal Calculated')
        disp(num2str(CPCT_calc_diag))
    end

    
    ERROR_CALC.CPC_calc_var = CPCT_calc;
    ERROR_CALC.CPC_calc_diag_var = CPCT_calc_diag;
    ERROR_CALC.CPC_infty_var = CPCT_exp;
    ERROR_CALC.CPC_infty_diag_var = CPCT_exp_diag;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
