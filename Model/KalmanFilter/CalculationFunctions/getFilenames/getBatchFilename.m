function [filename] = getBatchFilename(SIM,FLAG)
    if FLAG.InputMode == 5
        just_file = ['PRBS_Amp' num2str(FLAG.PRBSAmp) '_SwitchingTime' num2str(FLAG.Tswitch) '_SOC' num2str(FLAG.SOC) '_Ts' num2str(FLAG.Ts) '.mat'];
    else
    end

    filename = [FLAG.folderpathOptiBatch filesep just_file];
end