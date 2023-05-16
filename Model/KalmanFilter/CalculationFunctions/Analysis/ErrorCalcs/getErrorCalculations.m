function [ERROR_CALC,RESULTS] = getErrorCalculations(SIM,FLAG,N,P,RESULTS,ESTIMATOR)
N_samples = length(RESULTS.EST.VAR.t_soln);
idx = floor(FLAG.FractionOfData*N_samples);

fig = gcf;
NumFig = fig.Number;
if NumFig == 1
    close(1)
    NumFig = 0;
end


%% From DARE
for OO = 1:length(RESULTS.EST.PLANT.z_soln_ALL)
    if FLAG.EstimatorModel == 1 || FLAG.EST.SepHK == 0
        CPCT_DARE{OO}    = ESTIMATOR.C_des_ROM{1} * ESTIMATOR.P_infty{1} * ESTIMATOR.C_des_ROM{1}';
        CPCT_DARE_allROM = ESTIMATOR.C_des_ROM{1} * ESTIMATOR.P_infty{1} * ESTIMATOR.C_des_ROM{1}';
    else
        CPCT_DARE{OO} = ESTIMATOR.C_des_ROM{OO} * ESTIMATOR.P_infty{OO} * ESTIMATOR.C_des_ROM{OO}';
        CPCT_DARE_allROM(OO,OO) = CPCT_DARE{OO}(2,2);
    end
end


%% Calculate the error in the outputs (ASY)
if isfield(RESULTS.EST.ASY,'z_soln_ALL')
    for OO = 1:length(RESULTS.EST.PLANT.z_soln_ALL)
        plant_outputs = RESULTS.EST.PLANT.z_soln_ALL{OO};
        est_outputs   = RESULTS.EST.ASY.z_soln_ALL{OO};

        if FLAG.PlotError
            error_outputs = abs(est_outputs - plant_outputs);
            if FLAG.EST.SepHK
                figure(NumFig + OO)
                hold on
                plot(RESULTS.EST.VAR.t_soln   , error_outputs(2,:),'r','Linewidth',2,'DisplayName','Asymptotic')
                xline(idx,'Linewidth',2)
                xlabel('Time [s]')
                ylabel([RESULTS.Labels.unit{OO}])
                title([RESULTS.Labels.title{OO}  ' Output Error'])
                lgn = legend;
            else
                for i = 1:N.DesOut
                    figure(NumFig + i)
                    hold on
                    plot(RESULTS.EST.VAR.t_soln   , error_outputs(i,:),'r','Linewidth',2,'DisplayName','Asymptotic')
                    xline(idx,'Linewidth',2)
                    xlabel('Time [s]')
                    ylabel([RESULTS.Labels.unit{i}])
                    title([RESULTS.Labels.title{i}  ' Output Error'])
                    lgn = legend;
                end
            end
            
        end

        [N_outputs, ~] = size(plant_outputs);
        error_outputs  = est_outputs(:,idx:end)   - plant_outputs(:,idx:end);

        covar_outputs = nan(N_outputs);
        for i = 1:N_outputs
            for j = 1:N_outputs
                [covar_outputs(i,j)] = calcPopCovar(error_outputs(i,:), error_outputs(j,:));
            end
        end

        CPCT_calc_asy{OO} = covar_outputs;

        if FLAG.EstimatorModel == 1 || FLAG.EST.SepHK == 0
            CPCT_calc_asy_allROM = covar_outputs;
        else
            CPCT_calc_asy_allROM(OO,OO) = covar_outputs(2,2);
        end
    end
end


%% Calculate the error in the outputs (VAR)
if isfield(RESULTS.EST.VAR,'z_soln_ALL')
    for OO = 1:length(RESULTS.EST.PLANT.z_soln_ALL)
        plant_outputs = RESULTS.EST.PLANT.z_soln_ALL{OO};
        est_outputs   = RESULTS.EST.VAR.z_soln_ALL{OO};

        if FLAG.PlotError
            error_outputs = abs(est_outputs - plant_outputs);
            if FLAG.EST.SepHK
                figure(NumFig + OO)
                hold on
                plot(RESULTS.EST.VAR.t_soln   , error_outputs(2,:),'g','Linewidth',2,'DisplayName','Variable')
                xline(idx,'Linewidth',2)
                xlabel('Time [s]')
                ylabel([RESULTS.Labels.unit{OO}])
                title([RESULTS.Labels.title{OO}  ' Output Error'])
                lgn = legend;
            else
                for i = 1:N.DesOut
                    figure(NumFig + i)
                    hold on
                    plot(RESULTS.EST.VAR.t_soln   , error_outputs(i,:),'g','Linewidth',2,'DisplayName','Variable')
                    xline(idx,'Linewidth',2)
                    xlabel('Time [s]')
                    ylabel([RESULTS.Labels.unit{i}])
                    title([RESULTS.Labels.title{i}  ' Output Error'])
                    lgn = legend;
                end
            end
        end

        [N_outputs, ~] = size(plant_outputs);
        error_outputs = est_outputs(:,idx:end)   - plant_outputs(:,idx:end);

        covar_outputs = nan(N_outputs);
        for i = 1:N_outputs
            for j = 1:N_outputs
                [covar_outputs(i,j)] = calcPopCovar(error_outputs(i,:), error_outputs(j,:));
            end
        end

        CPCT_calc_var{OO} = covar_outputs;

        if FLAG.EstimatorModel == 1 || FLAG.EST.SepHK == 0
            CPCT_calc_var_allROM = covar_outputs;
        else
            CPCT_calc_var_allROM(OO,OO) = covar_outputs(2,2);
        end
    end
end


%% Some Calcs
CPCT_DARE_allROM_diag  = diag(CPCT_DARE_allROM)';
if exist('CPCT_calc_asy_allROM','var')
    CPCT_calc_asy_allROM_diag = diag(CPCT_calc_asy_allROM)';
end
if exist('CPCT_calc_var_allROM','var')
    CPCT_calc_var_allROM_diag = diag(CPCT_calc_var_allROM)';
end


%% Display Results to Command Window
    if FLAG.Analysis.dispResults
        if 1 && FLAG.DoAsy
            %% Asymptotic
            disp(newline)
            disp('____________________')
            disp('Asymptotic Estimator')
            disp('____________________')
            
            % Each Output
            if 1
                disp(newline)
                disp('Each Output')
                disp('-----------')
                disp(newline)
                for OO = 1:length(CPCT_calc_asy)
                    disp(['Output ' RESULTS.Labels.title{OO}])
                    %disp(newline)
                    disp('   DARE')
                    disp(num2str(CPCT_DARE{OO}))
                    disp('   Calculated')
                    disp(num2str(CPCT_calc_asy{OO}))
                    disp(newline)
                end
            end
    
            % Combined Outputs
            if 1
                disp(newline)
                disp('Combined Output')
                disp('---------------')
                disp('   DARE')
                disp(num2str(CPCT_DARE_allROM))
                disp(newline)
                disp('   Calculated')
                disp(num2str(CPCT_calc_asy_allROM))
                disp(newline)
                disp('With difference of')
                disp(num2str(CPCT_calc_asy_allROM - CPCT_DARE_allROM))
                disp(newline)
                disp('Ratio of error')
                disp(num2str(abs((CPCT_calc_asy_allROM - CPCT_DARE_allROM)./CPCT_DARE_allROM)))
                disp(newline)
            end
    
            % Just the diagonal of all CPCT
            if 1
                disp(newline)
                disp('Diagonals')
                disp('---------')
                disp('CPC^T diagonal DARE')
                disp(num2str(CPCT_DARE_allROM_diag))
                disp(newline)
                disp('CPC^T diagonal Calculated')
                disp(num2str(CPCT_calc_asy_allROM_diag))
            end
        end

    %% Variable Estimator
        if 1 && FLAG.DoVar
            disp(newline)
            disp('__________________')
            disp('Variable Estimator')
            disp('__________________')
            
            % Each Output
            if 1
                disp(newline)
                disp('Each Output')
                disp('-----------')
                disp(newline)
                for OO = 1:length(CPCT_calc_var)
                    disp(['Output ' RESULTS.Labels.title{OO}])
                    %disp(newline)
                    disp('   DARE')
                    disp(num2str(CPCT_DARE{OO}))
                    disp('   Calculated')
                    disp(num2str(CPCT_calc_var{OO}))
                    disp(newline)
                end
            end
    
            % Combined Outputs
            if 1
                disp(newline)
                disp('Combined Output')
                disp('---------------')
                disp('   DARE')
                disp(num2str(CPCT_DARE_allROM))
                disp(newline)
                disp('   Calculated')
                disp(num2str(CPCT_calc_var_allROM))
                disp(newline)
                disp('With difference of')
                disp(num2str(CPCT_calc_var_allROM - CPCT_DARE_allROM))
                disp(newline)
                disp('Ratio of error')
                disp(num2str(abs((CPCT_calc_var_allROM - CPCT_DARE_allROM)./CPCT_DARE_allROM)))
                disp(newline)
            end
    
            % Just the diagonal of all CPCT
            if 1
                disp(newline)
                disp('Diagonals')
                disp('---------')
                disp('CPC^T diagonal DARE')
                disp(num2str(CPCT_DARE_allROM_diag))
                disp(newline)
                disp('CPC^T diagonal Calculated')
                disp(num2str(CPCT_calc_var_allROM_diag))
            end
    
        end
    end

%% Save Outputs
    ERROR_CALC.CPC_DARE_all      = CPCT_DARE_allROM;
    ERROR_CALC.CPC_DARE_all_diag = CPCT_DARE_allROM_diag;

    if FLAG.DoAsy
        ERROR_CALC.CPC_calc_all_asy      = CPCT_calc_asy_allROM;
        ERROR_CALC.CPC_calc_all_diag_asy = CPCT_calc_asy_allROM_diag;
    end
    
    if FLAG.DoVar
        ERROR_CALC.CPC_calc_all_var      = CPCT_calc_var_allROM;
        ERROR_CALC.CPC_calc_all_diag_var = CPCT_calc_var_allROM_diag;
    end


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
