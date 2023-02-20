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
    if ~exist(ESTIMATOR)
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
    Kalman_plotFnc(RESULTS,N)
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%% OOOOOOOOOOOLD %%%%%%%%%%%%%%%%%%%%%%%%%% 
% %% Check if File Exists
% if isfile(SIM.save_filename)
%     if SIM.Analysis.OverwriteData
%         disp('Simulation already exists, deleting data')
%         delete(SIM.save_filename)
%     else
%         disp('Simulation already exists, loading data')
%         load(SIM.save_filename)
%     end
% end


%% Initialization
% if ~ SIM.AnalysisComplete.Initialization
%     [SIM,FLAG] = inputs(SIM,FLAG);

% %     save(SIM.save_filename,"FLAG","P","N","SIM")
%     
% %     SIM.AnalysisComplete.Initialization = 1;
% end

%% No Noise Comparison
%     old_C_mode = FLAG.C_mode;
%     FLAG.C_mode = 5;

    
%     FLAG.C_mode = old_C_mode;
    
%     save(SIM.save_filename,"FLAG","P","N","SIM","RESULTS")

%%
 %&& ~SIM.AnalysisComplete.NoisyPlant


%     SIM.AnalysisComplete.NoisyPlant = 1;
%     
%     % Save Overall
%         save(SIM.save_filename,"FLAG","P","N","SIM","RESULTS")
% 
    % Save Just Plant Data

%% Estimator
%&&  ~SIM.AnalysisComplete.Estimator


%     SIM.AnalysisComplete.Estimator = 1;
%     save(SIM.save_filename,"FLAG","P","N","SIM","RESULTS","ESTIMATOR")

%% Error Calculations
 %&& ~ SIM.AnalysisComplete.Est_Error_calc


    
%     SIM.AnalysisComplete.Est_Error_calc = 1;

%     save(SIM.save_filename,"FLAG","P","N","SIM","RESULTS","ESTIMATOR","ERROR_CALC")

%% Generate Comparison Data
%&& ~SIM.AnalysisComplete.GenComparData


%     SIM.AnalysisComplete.GenComparData = 1;

%     save(SIM.save_filename,"FLAG","P","N","SIM","RESULTS","ESTIMATOR","ERROR_CALC")


%% Compare SVD to P_infty
%&& ~SIM.AnalysisComplete.ComparSVD2Pinf


%     
% %     SIM.AnalysisComplete.ComparSVD2Pinf = 1;
% 
%     save(SIM.save_filename,"FLAG","P","N","SIM","RESULTS","ESTIMATOR","ERROR_CALC")
