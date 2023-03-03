function filename = getROMErrorFilename(FLAG)
    just_file = ['ROMError_SOC' num2str(FLAG.SOC) '_Ts' num2str(FLAG.Ts) '.mat'];
    filename = [FLAG.folderpathROMError filesep just_file];
end