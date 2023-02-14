
function [overall_filename] = getSlinkfilename(FLAG,SIM)
    switch FLAG.SlinkModel
        case 3 % 3 state plant with state noise
            type_str = 'NPlant3';
        case 5
            type_str = 'NPlant5';
    end
    
    switch FLAG.InputMode
        case 1
            input_str = '_Step';
        case 2
            input_str = '_Sine';
            fq_str    = [num2str(SIM.fq) 'Hz'];
            input_str = [input_str fq_str];
    end

    Q_str = num2str(SIM.Q_0);
    R_str = num2str(SIM.R_0);
    Ts_str = num2str(SIM.Ts);
    switch FLAG.C_mode
        case 1
            C_str = num2str(FLAG.C_mode);
        case 2
            C_str = num2str(FLAG.C_mode);
        case 3
            C_str = num2str(FLAG.C_mode);
        case 4
            C_str = num2str('All');
    end
    
    
    Folder = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanTesting\DAE_Circuit\SoftwareStyle\Results\PlantData';
    filename = [type_str input_str '_Q' Q_str '_R' R_str '_Ts' Ts_str '_C' C_str '.mat'];
    overall_filename = [Folder filesep filename];

end