%% get Simulink Filename
function [overall_filename] = getSlinkfilename(FLAG,SIM)
   
    switch FLAG.InputMode
        case 1
            input_str = 'Step';
        case 2
            input_str = 'Sine';
            fq_str    = [num2str(SIM.fq) 'Hz'];
            input_str = [input_str fq_str];
    end

    Q_str  = num2str(SIM.Q_0);
    R_str  = num2str(SIM.R_0);
    Ts_str = num2str(SIM.Ts);
%     switch FLAG.C_mode
%         case 5
%             C_str = num2str('All');
%         otherwise
%             C_str = num2str(FLAG.C_mode);
%     end
    SOC_str = num2str(FLAG.SOC);
    
    
    Folder = [pwd filesep 'Results' filesep 'PlantData'];
    filename = [input_str '_Q' Q_str '_R' R_str '_Ts' Ts_str '_SOC' SOC_str '.mat'];
    overall_filename = [Folder filesep filename];

end