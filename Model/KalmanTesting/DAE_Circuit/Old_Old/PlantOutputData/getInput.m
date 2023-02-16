%% Get Ramp Voltage Input
function [V_in] = getInput(t,SIM)

    if t <= SIM.time_rest
        V_in = SIM.V_init;
    elseif t > SIM.time_rest && t < (SIM.time_rest + SIM.time_ramp)
        V_in = (SIM.V_step - SIM.V_init)/SIM.time_ramp*(t-SIM.time_rest) + SIM.V_init;
    else
        V_in = SIM.V_step;
    end

end