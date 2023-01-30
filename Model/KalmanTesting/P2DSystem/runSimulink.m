%% Run P2D Simulink
%
%
% Slink will hold things such as tfinal, (A,B,C)_ROM, (A,B,C)_est, x0,
% x0_est C, Ts, Iden, P_k_0
function [out] = runSimulink(InputSignal,Noise,Slink,LC_param,sim_in,localFLAG)
%     t_final = InputSignal(end,1);
    t_final = InputSignal.t_final;
    if localFLAG.InputType == 1
        step_time = InputSignal.step_time;
        step_amp = InputSignal.step_amp;
    else
        sine_amp  = InputSignal.sine_amp;
        sine_freq = InputSignal.omega;
    end

    Q = Noise.Q;
    R = Noise.R;
    Ts_noise = Noise.Ts_noise;
    pow_Q = Noise.pow_Q;
    pow_R = Noise.pow_R;
  
    C_plant  = Slink.C_plant;
    A_r      = Slink.A_r;
    B_r      = Slink.B_r;
    C_r_meas = Slink.C_r_meas;
    x0_ROM   = Slink.x0_ROM;
    x0_est   = Slink.x0_est;
    P_k_0    = Slink.P_k_0;
    Iden     = Slink.Iden;
    K_infty  = Slink.K_infty;

    
    M_cross = LC_param.M_cross;
    U_NV_NT = LC_param.U_NV_NT;
    K_plant = LC_param.K;
    
    AN    = sim_in.AN;
    CA    = sim_in.CA;
    SEP   = sim_in.SEP;
    EL    = sim_in.EL;
    SIM   = sim_in.SIM;
    CONS  = sim_in.CONS;
    P     = sim_in.P;
    N     = sim_in.N;
    FLAG  = sim_in.FLAG;
    PROPS = sim_in.PROPS;
    SV_IC = sim_in.SIM.SV_IC;

    Ts = sim_in.SIM.Ts;

%% Function Handle Conversion
    AN.EqPotentialHandle = func2str(AN.EqPotentialHandle);
    AN.i_oHandle         = func2str(AN.i_oHandle);
    AN.sigmaHandle       = func2str(AN.sigmaHandle);
    AN.D_oHandle         = func2str(AN.D_oHandle);
    
    CA.EqPotentialHandle = func2str(CA.EqPotentialHandle);
    CA.i_oHandle         = func2str(CA.i_oHandle);
    CA.sigmaHandle       = func2str(CA.sigmaHandle);
    CA.D_oHandle         = func2str(CA.D_oHandle);
    
    EL.tf_numHandle      = func2str(EL.tf_numHandle);
    EL.ActivityHandle    = func2str(EL.ActivityHandle);
    EL.D_o_Li_ionHandle  = func2str(EL.D_o_Li_ionHandle);
    EL.kappaHandle       = func2str(EL.kappaHandle);

    SIM.ControllerHandle = func2str(SIM.ControllerHandle);
    
    SIM = rmfield(SIM,'fsolve_options');
    SIM = rmfield(SIM,'Controller_MO_File');

% Model Parameters
    if localFLAG.InputType == 1
        mdl = 'P2D_Overall';
    else
        mdl = 'P2D_Overall_sine';
    end
    load_system(mdl)
    in = Simulink.SimulationInput(mdl);
    in = in.setModelParameter('StartTime','0','StopTime',num2str(t_final));
    
    mdlWks = get_param(in,'ModelWorkspace');
    if localFLAG.InputType == 1
        assignin(mdlWks,'step_time' ,step_time)
        assignin(mdlWks,'step_amp' ,step_amp)
    else
        assignin(mdlWks,'sine_freq' ,sine_freq)
        assignin(mdlWks,'sine_amp' ,sine_amp)
    end
    assignin(mdlWks,'Q'     ,Q)
    assignin(mdlWks,'R'     ,R)
    assignin(mdlWks,'Ts'     ,Ts)
    assignin(mdlWks,'Ts_noise'     ,Ts_noise)
    assignin(mdlWks,'pow_Q'     ,pow_Q)
    assignin(mdlWks,'pow_R'     ,pow_R)
    assignin(mdlWks,'C_plant',C_plant)
    assignin(mdlWks,'M_cross',M_cross)
    assignin(mdlWks,'U_NV_NT',U_NV_NT)
    assignin(mdlWks,'K_plant',K_plant)
    assignin(mdlWks,'AN',AN)
    assignin(mdlWks,'CA',CA)
    assignin(mdlWks,'SEP',SEP)
    assignin(mdlWks,'EL',EL)
    assignin(mdlWks,'SIM',SIM)
    assignin(mdlWks,'CONS',CONS)
    assignin(mdlWks,'P',P)
    assignin(mdlWks,'N',N)
    assignin(mdlWks,'FLAG',FLAG)
    assignin(mdlWks,'PROPS',PROPS)
    assignin(mdlWks,'SV_IC',SV_IC)
    
    assignin(mdlWks,'A_ROM' ,A_r)
    assignin(mdlWks,'B_ROM' ,B_r)
    assignin(mdlWks,'C_ROM' ,C_r_meas)
    assignin(mdlWks,'x0_ROM',x0_ROM)

    assignin(mdlWks,'A_est' ,A_r)
    assignin(mdlWks,'B_est' ,B_r)
    assignin(mdlWks,'C_est' ,C_r_meas)
    assignin(mdlWks,'x0_est',x0_est)

    assignin(mdlWks,'P_k_0',P_k_0)
    assignin(mdlWks,'Iden' ,Iden)
    assignin(mdlWks,'K_infty',K_infty)

    out = sim(in);
%     save_system
%     close_system

end