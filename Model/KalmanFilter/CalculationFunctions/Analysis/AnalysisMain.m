%% Analysis Main
% Determines which of the analysis needs to be performed

function AnalysisMain(FLAG)
%% Initialization
    [SIM,N,P,FLAG,RESULTS] = init(FLAG);


%% Perform Analysis
%% No Noise Comparison
if FLAG.Analysis.NoNoiseCompare
    [SIM,FLAG,N,P,RESULTS] = NoNoiseComparisonMain(SIM,FLAG,N,P,RESULTS);
end


%% get DARE
if FLAG.Analysis.getDARE
    [RESULTS] = getDARE(SIM,FLAG,N,P,RESULTS);
end


%% Noisy Plant
if FLAG.Analysis.NoisyPlant
    [RESULTS] = RunSimulink(SIM,N,P,FLAG,RESULTS);

    % Save Data from Simulink
        Plant_Data = RESULTS.Slink_plant;
        [Slink_filename] = getSlinkfilename(FLAG,SIM);
        save(Slink_filename,"FLAG","P","N","SIM","Plant_Data")
end


%% Estimator
if FLAG.Analysis.Estimator 
    [plant_filename] = getSlinkfilename(FLAG,SIM);
    FLAG.Analysis.PlotImp = 0;
    [ESTIMATOR,RESULTS] = PerformEstimation(plant_filename,SIM,FLAG,N,P,RESULTS);

    %% Save Data for plots
        % save('F:\TylerFiles\~PhDWork\Presentations\Seminar\Spring2023\Pictures\EstimationResults\EstimationResultsData.mat','RESULTS')

    %% Plot Kalman Gain
    if FLAG.PLOT.K_k_gain
        [~,~,N_steps] = size(ESTIMATOR.K_k{1});
        norm_vec = nan(1,N_steps);
        for i = 1:N_steps % get norm
           norm_vec(i) = norm(ESTIMATOR.K_k{1}(:,:,i));
        end

        figure
        plot(RESULTS.EST.PLANT.t_soln , norm_vec, 'k','Linewidth',2)
        title('Norm of Kalman Gain')
        xlabel('Time [s]')
        ylabel('|K_k|')

        min_t = RESULTS.EST.PLANT.t_soln(end-50);
        max_t = RESULTS.EST.PLANT.t_soln(end);
        figure
        plot(RESULTS.EST.PLANT.t_soln , norm_vec, 'k','Linewidth',2)
        title('Norm of Kalman Gain: Last 50 Points')
        xlabel('Time [s]')
        ylabel('|K_k|')
        xlim([min_t,max_t])
    end

end


%% Error Calculations
if FLAG.Analysis.Est_Error_calc
    if ~exist('ESTIMATOR')
        ESTIMATOR = struct();
    end
    [ERROR_CALC,RESULTS] = getErrorCalculations(SIM,FLAG,N,P,RESULTS,ESTIMATOR);
end


%% Generate Comparison Data
if FLAG.Analysis.GenComparData 
    generateComparisonData(SIM,FLAG,N,P,RESULTS) % ,ESTIMATOR,ERROR_CALC
end


%% Compare SVD to P_infty
if FLAG.Analysis.ComparSVD2Pinf 
    CompareData(FLAG);
end


%% Plot Results
if FLAG.PLOT.PlotResults
    Kalman_plotFnc(RESULTS,N,SIM,FLAG)
end


end

