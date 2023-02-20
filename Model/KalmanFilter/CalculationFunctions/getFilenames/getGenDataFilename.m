function filename = getGenDataFilename(FLAG)
    filename = ['CompData_Q' num2str(FLAG.Q_0) '_R' num2str(FLAG.R_0) '_Ts' num2str(FLAG.Ts) '_SOC' num2str(FLAG.SOC) '.mat'];
    if FLAG.QMode == 1 % Input Q
        path_str = [pwd filesep 'Results' filesep 'ComparisonData' filesep 'InputQ'];
    else % State Q
        path_str = [pwd filesep 'Results' filesep 'ComparisonData' filesep 'StateQ'];
    end
    filename = [path_str filesep filename];
end