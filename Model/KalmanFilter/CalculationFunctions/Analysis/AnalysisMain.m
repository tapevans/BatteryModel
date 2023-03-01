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
    [ESTIMATOR,RESULTS] = PerformEstimation(plant_filename,SIM,FLAG,N,P,RESULTS);
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
    Kalman_plotFnc(RESULTS,N,SIM)
end


end

