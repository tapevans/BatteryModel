%% Initialization File
function [SIM,N,P,FLAG] = init(SIM,FLAG)
%% Pointers
    i = 1;
    P.V_1   = i; i = i+1;
    P.V_2   = i; i = i+1;
    P.V_3   = i; i = i+1;
    P.V_s   = i; i = i+1;
    P.i_PS  = i; i = i+1;


%% N Variables
%     if FLAG.NStates == 3
%         [~,C_m]=getSS_C3(SIM,P,FLAG);
%     else
%         [~,C_m]=getSS_C5(SIM,P,FLAG);
%     end
%     [N.measur , ~]  = size(C_m);
switch FLAG.C_mode
    case 1
        N.measur = 1;
    case 2
        N.measur = 1;
    case 3
        N.measur = 1;
    case 4
        N.measur = 3;
end
    N.inputs = 1;


%% Mass Matrix
% 3 State System
    Mass = zeros(3);
    
    Mass(P.V_1 , P.V_1 ) =  SIM.C_1;
    Mass(P.V_3 , P.V_3 ) =  SIM.C_4;

    SIM.Mass3 = Mass;

    SIM.algb_idx3 = [P.V_2];
    SIM.diff_idx3 = [P.V_1 , P.V_3];

% 5 State System
    Mass = zeros(5);
     
    Mass(P.V_1 , P.V_1) =  SIM.C_1;
    Mass(P.V_1 , P.V_s) = -SIM.C_1;
    Mass(P.V_s , P.V_1) = -SIM.C_1;
    Mass(P.V_s , P.V_s) =  SIM.C_1;
    Mass(P.V_3 , P.V_3) =  SIM.C_4;

    SIM.Mass5 = Mass;

    SIM.algb_idx5 = [P.V_2 , P.i_PS];
    SIM.diff_idx5 = [P.V_1 , P.V_3, P.V_s];


%% Input Signal Parameters
    SIM.V_step = SIM.V_init + SIM.V_step_height;
    
%     SIM.time_rest = 8*SIM.tc+SIM.Ts;
    SIM.time_rest = 0;
    SIM.time_ramp = 2*SIM.tc;
    
    SIM.t_final_sim = SIM.t_final + 2*SIM.Ts;
    SIM.t_vec = 0:SIM.tc:SIM.t_final_sim;


%% Initial Conditions
% Initial Steady State
    R_tot = SIM.R_1 + SIM.R_2 + SIM.R_3 + SIM.R_4;
    i_PS  = SIM.V_init / R_tot;
    V_1   = SIM.V_init - SIM.R_1 * i_PS;
    V_2   = V_1        - SIM.R_2 * i_PS;
    V_3   = V_2        - SIM.R_3 * i_PS;
    delV_1= SIM.V_init - V_1;
    delV_2= V_1        - V_2;
    delV_3= V_2        - V_3;
    delV_4= V_3        - 0;

% Initial Step
    V_eff = SIM.V_step - delV_1 -delV_4;
    R_eff = SIM.R_2 + SIM.R_3;
    i_eff = V_eff/R_eff;
    delV_2_eff = i_eff*SIM.R_2;

% 3 State System
    SIM.x_0_3        = zeros(P.V_3,1);
    SIM.x_0_3(P.V_1) = SIM.V_step - delV_1;
    SIM.x_0_3(P.V_3) = delV_4;
    SIM.x_0_3(P.V_2) = (SIM.V_step - delV_1) - delV_2_eff ;

% 5 State System
    SIM.x_0_5 = zeros(5,1);
    
    SIM.x_0_5(P.V_s)  = SIM.V_init;
    SIM.x_0_5(P.i_PS) = i_PS;
    SIM.x_0_5(P.V_1)  = V_1;
    SIM.x_0_5(P.V_2)  = V_2;
    SIM.x_0_5(P.V_3)  = V_3;


%% Calculate Noise Power
    SIM.Q3i = SIM.Q_0 * eye(N.inputs);
    SIM.Q5i = SIM.Q_0 * eye(N.inputs);
    
    Q       = SIM.Q_0 * eye(3);
    SIM.Q3s = (diag(Q));
    Q       = SIM.Q_0 * eye(5);
    SIM.Q5s = (diag(Q))';
    
    R     = SIM.R_0 * eye(N.measur);
    SIM.R = (diag(R));
    
    SIM.pow_Q3i = SIM.tc * SIM.Q3i;
    SIM.pow_Q5i = SIM.tc * SIM.Q5i;
    SIM.pow_Q3s = SIM.tc * SIM.Q3s;
    SIM.pow_Q5s = SIM.tc * SIM.Q5s;
    SIM.pow_R   = SIM.tc * SIM.R;

end