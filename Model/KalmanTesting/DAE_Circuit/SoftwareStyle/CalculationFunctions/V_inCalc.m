%% V_inCalc
function [V_in] = V_inCalc(t_vec,SIM,FLAG)

switch FLAG.InputType
    case 1 % Ramp Step
        for i = 1:length(t_vec)
            t = t_vec(i);
            if t <= SIM.time_rest
                V_in(i) = SIM.V_init;
            elseif t > SIM.time_rest && t < (SIM.time_rest + SIM.time_ramp)
                V_in(i) = (SIM.V_step - SIM.V_init)/SIM.time_ramp*(t-SIM.time_rest) + SIM.V_init;
            else
                V_in(i) = SIM.V_step;
            end

        end

    case 2 % Sine

    case 3 % Step
        V_in = SIM.V_step*ones(size(t_vec));
end

end


