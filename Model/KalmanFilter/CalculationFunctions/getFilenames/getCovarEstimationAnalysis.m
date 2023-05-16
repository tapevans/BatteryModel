function filename = getCovarEstimationAnalysis(FLAG , QQ , RR)
    just_file = ['CovarEst_Combined_SOC' num2str(FLAG.SOC) '_Ts' num2str(FLAG.Ts) , '_Q' , num2str(QQ) , '_R' , num2str(RR) , '.mat'];
    filename = [FLAG.folderpathCovarEst filesep just_file];
end