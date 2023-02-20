function [sys_CT , sys_DT] = getSS_System(SIM,N,P,FLAG)
    % Get Filename
    filename = getSSFilename(FLAG);
    EIS_data = load(filename);
    
    % Produce SS
    D = zeros(N.DesOut,N.inputs);
    sys = dss(EIS_data.A , EIS_data.B , SIM.OutputMatrix , D , EIS_data.SIM.M);
    sys_CT = ss(sys,'explicit');
    sys_DT = c2d(sys_CT,SIM.Ts);
end