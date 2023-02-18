%% Analysis Main
% Determines which of the analysis needs to be performed

function AnalysisMain(SIM,FLAG)
%% Check if File Exists
if isfile(SIM.save_filename)
    if SIM.Analysis.OverwriteData
        disp('Simulation already exists, deleting data')
        delete(SIM.save_filename)
    else
        disp('Simulation already exists, loading data')
        load(SIM.save_filename)
    end
end

%% Perform Analysis
%% Initialization
if ~ SIM.AnalysisComplete.Initialization
    [SIM,FLAG] = inputs(SIM,FLAG);
    [SIM,N,P,FLAG,RESULTS] = init(SIM,FLAG);
    save(SIM.save_filename,"FLAG","P","N","SIM")
    
    SIM.AnalysisComplete.Initialization = 1;
end


%% No Noise Comparison
if FLAG.Analysis.NoNoiseCompare
    old_C_mode = FLAG.C_mode;
    FLAG.C_mode = 5;

    [SIM,FLAG,N,P,RESULTS] = NoNoiseComparisonMain(SIM,FLAG,N,P,RESULTS);
    
    FLAG.C_mode = old_C_mode;
    
    save(SIM.save_filename,"FLAG","P","N","SIM","RESULTS")
end


%% Noisy Plant
if FLAG.Analysis.NoisyPlant && ~SIM.AnalysisComplete.NoisyPlant
    [RESULTS] = RunSimulink(SIM,N,P,FLAG,RESULTS);

    SIM.AnalysisComplete.NoisyPlant = 1;
    
    % Save Overall
        save(SIM.save_filename,"FLAG","P","N","SIM","RESULTS")

    % Save Just Plant Data
        Plant_Data = RESULTS.Slink_plant;
        [overall_filename] = getSlinkfilename(FLAG,SIM);
        save(overall_filename,"FLAG","P","N","SIM","Plant_Data")
end


%% Estimator
if FLAG.Analysis.Estimator && ~SIM.AnalysisComplete.Estimator
    [plant_filename] = getSlinkfilename(FLAG,SIM);
    [ESTIMATOR,RESULTS] = PerformEstimation(plant_filename,SIM,FLAG,N,P,RESULTS);

    SIM.AnalysisComplete.Estimator = 1;
    save(SIM.save_filename,"FLAG","P","N","SIM","RESULTS","ESTIMATOR")
end


%% Error Calculations
if FLAG.Analysis.Est_Error_calc && ~ SIM.AnalysisComplete.Est_Error_calc
    [ERROR_CALC,RESULTS] = getErrorCalculations(SIM,FLAG,N,P,RESULTS,ESTIMATOR);
    
    SIM.AnalysisComplete.Est_Error_calc = 1;

    save(SIM.save_filename,"FLAG","P","N","SIM","RESULTS","ESTIMATOR","ERROR_CALC")
end


%% Generate Comparison Data
if FLAG.Analysis.GenComparData && ~SIM.AnalysisComplete.GenComparData
    generateComparisonData(SIM,FLAG,N,P,RESULTS,ESTIMATOR,ERROR_CALC)

    SIM.AnalysisComplete.GenComparData = 1;

    save(SIM.save_filename,"FLAG","P","N","SIM","RESULTS","ESTIMATOR","ERROR_CALC")
end


%% Compare SVD to P_infty
if FLAG.Analysis.ComparSVD2Pinf && ~SIM.AnalysisComplete.ComparSVD2Pinf
    CompareData(FLAG);
    
%     SIM.AnalysisComplete.ComparSVD2Pinf = 1;

    save(SIM.save_filename,"FLAG","P","N","SIM","RESULTS","ESTIMATOR","ERROR_CALC")
end


%% Plot Results
if FLAG.PLOT.PlotResults
    MSD_plotFnc(SIM.save_filename)
end


end

