function filename = getImpulseFilename(FLAG)
    Ts_vec_str = {'0.1',
    '0.121152765862859',
    '0.146779926762207',
    '0.177827941003892',
    '0.215443469003188',
    '0.261015721568254',
    '0.316227766016838',
    '0.383118684955729',
    '0.464158883361278',
    '0.562341325190349',
    '0.681292069057961',
    '0.825404185268018',
    '1.0',
    '1.21152765862859',
    '1.46779926762207',
    '1.77827941003892',
    '2.15443469003188',
    '2.61015721568254',
    '3.16227766016838',
    '3.83118684955729',
    '4.64158883361278',
    '5.62341325190349',
    '6.81292069057961',
    '8.25404185268018',
    '10.0'};

    for i = 1:length(Ts_vec_str)
        Ts_vec_str2num(i) = str2num(Ts_vec_str{i});
    end

    [~,idx] = min(abs(Ts_vec_str2num - FLAG.Ts));



    just_file = ['ObservabilityTest_KPCont_DTImpulseTs' num2str(Ts_vec_str{idx}) 'SOC' num2str(FLAG.SOC) '.mat'];
    filename = [FLAG.folderpath filesep just_file];
end