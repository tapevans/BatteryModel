%% get Simulink Filename
function [overall_filename] = getSlinkfilename(FLAG,SIM)
   
    switch FLAG.InputMode
        case 1
            input_str = 'Step';
        case 2
            input_str = 'Sine';
            fq_str    = [num2str(SIM.fq) 'Hz'];
            input_str = [input_str fq_str];
        case 5
            input_str   = 'PRBS';
            PRBSAmp_str = ['_Amp' num2str(FLAG.PRBSAmp)];
            Tswitch_str = ['_Tswitch' num2str(FLAG.Tswitch)];
            input_str   = [input_str PRBSAmp_str Tswitch_str];
    end

    Q_str  = num2str(SIM.Q_0);
    R_str  = num2str(SIM.R_0);
    Ts_str = num2str(SIM.Ts);
    SOC_str = num2str(FLAG.SOC);
    
    Folder = [pwd filesep 'Results' filesep 'PlantData'];
    filename = [input_str '_Q' Q_str '_R' R_str '_Ts' Ts_str '_SOC' SOC_str '.mat'];
    overall_filename = [Folder filesep filename];

end