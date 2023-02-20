%% Test Getting Noise for Slink
%     Q_0_vec = [1e-6];
%     R_0_vec = [1e-6];
%     Ts_vec  = [1];
%     SOC_vec = [50];
% 
%     FLAG.Q_0 = QQ; % Process Noise
%     FLAG.R_0 = RR; % Measurement Noise
%     FLAG.Ts = TT;  % Sampling Rate (When to save outputs)
%     FLAG.SOC = SS; % State of Charge
% 
%     SIM.Q_0 = FLAG.Q_0;
%     SIM.R_0 = FLAG.R_0;
%     SIM.Ts  = FLAG.Ts;
%     SIM.tc  = FLAG.Ts;

    SIM.Q_0 = [1e-6];
    SIM.R_0 = [1e-6];
    SIM.Ts  = [1];
    SIM.SOC = [50];

    SIM.tc = SIM.Ts;

    SIM.Qi = SIM.Q_0 * eye(N.inputs);
    
    Q       = SIM.Q_0 * eye(N.states);
    SIM.Qs = (diag(Q));
    
    R     = SIM.R_0 * eye(N.measur);
    SIM.R = (diag(R));
    
    SIM.pow_Qi = SIM.tc * SIM.Qi;
    SIM.pow_Qs = SIM.tc * SIM.Qs;
    SIM.pow_R  = SIM.tc * SIM.R;