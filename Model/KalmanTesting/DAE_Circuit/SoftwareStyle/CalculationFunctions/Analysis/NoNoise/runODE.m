%% runODE
function [t_soln, x_soln, z_soln] = runODE(SIM,N,P,FLAG)
%% InputSignalMode
%  1) From code (From V_inCalc)
%  2) From function input (For Simulink)
    FLAG.InputSignalMode = 1;

    switch FLAG.StateMode
        case 3
            x_0 = SIM.x_0_3;
            Mass = SIM.Mass3;
            [A,B,C_m,~] = getAll_SS3(SIM,N,P,FLAG);
        case 5
            x_0 = SIM.x_0_5;
            Mass = SIM.Mass5;
            [A,B,C_m,~] = getAll_SS5(SIM,N,P,FLAG);
    end

    if FLAG.InputMode == 4
        tspan = SIM.t_vec_imp;
    else
        tspan = [0,SIM.t_final_sim];
    end
    Tol.Rel = 1e-4;
    Tol.Abs = 1e-7;
    options_DAE = odeset('RelTol' ,Tol.Rel,      ...
                         'AbsTol' ,Tol.Abs,      ...
                         'Mass'   ,Mass);
    V_in = SIM.V_init;
    
    switch FLAG.StateMode
        case 3
            [t_soln,x_soln] = ode15s(@(t,SV) GovnEqn_3(t,SV,SIM,N,P,FLAG,A,B,V_in),tspan,x_0,options_DAE);
        case 5
            [t_soln,x_soln] = ode15s(@(t,SV) GovnEqn_5(t,SV,SIM,N,P,FLAG,A,B,V_in),tspan,x_0,options_DAE);
    end
    z_soln = (C_m*x_soln')';


end