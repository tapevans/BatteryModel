function filename = getODEStepFilename(FLAG)
    just_file = ['ObservabilityTest_KPCont_Relax_StepSOC' num2str(FLAG.SOC) '.mat'];
    filename = [FLAG.folderpath filesep just_file];
end