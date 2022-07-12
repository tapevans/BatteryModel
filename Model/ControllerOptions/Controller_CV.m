%% Controller_CV
% Modes of Operation (MO)
% 0) Simulation Ended
% 1) Constant Current
% 2) Constant Voltage
% 3) Relaxation

function [i_user] = Controller_CV(SV, P , SIM)
% Parameters
% 	K_p = 1200;  % Too Aggressive
%     K_p = 100; % Too Slow
%     K_p = 700; % Too Aggressive
%     K_p = 500; % Too Aggressive
%     K_p = 400; % Aggressive
    K_p = 275; % Too Slow
%     K_p = 300; % Slightly unstable still
    
	K_i = 0;
	K_d = 0;
    C_rate_min = 1/20;
    C_rate_max = 2;
    
%% Cell Voltage Calc
	cell_voltage = SV(P.phi_ed,end) - SV(P.phi_ed,1);

%% Error Calculation
	voltage_ref = SIM.Controller_MO_File(SIM.current_MO_step).Volt_ref;
	error = voltage_ref - cell_voltage;

%% Output Current Calculation
    i_calc = K_p * -error;
%            + K_i * trapz(time_history,error_history) ...
%            + K_d * (error_history(end) - error_history(end -1))/(time_history(end) - time_history(end -1));            

%% Saturation Limits
    i_min = C_rate_min * SIM.Cell_Cap / SIM.A_c;
    i_max = C_rate_max * SIM.Cell_Cap / SIM.A_c;
    if abs(i_calc) < i_min % if i_calc is smaller than min limit, set to zero
        i_user = 0;
    elseif abs(i_calc) > i_max % if i_calc is larger than max limit, saturate to i_max
        if i_calc < 0
            i_user = -i_max;
        else
            i_user =  i_max;
        end
    else
        i_user = i_calc;
    end

end