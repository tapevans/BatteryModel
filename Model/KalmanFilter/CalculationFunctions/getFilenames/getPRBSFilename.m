function filename = getPRBSFilename(FLAG)
    just_file = ['PRBS_Sims_PRBS_Amp' num2str(FLAG.PRBSAmp) '_SOC' num2str(FLAG.SOC) '_SwitchingTime' num2str(FLAG.Tswitch) '.mat'];
    filename = [FLAG.folderpathPRBS filesep just_file];
end