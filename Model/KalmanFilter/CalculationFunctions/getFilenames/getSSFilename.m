function filename = getSSFilename(FLAG)
    just_file = ['ObservabilityTest_SS_EIS_SOC' num2str(FLAG.SOC) '.mat'];
    filename = [FLAG.folderpath filesep just_file];
end