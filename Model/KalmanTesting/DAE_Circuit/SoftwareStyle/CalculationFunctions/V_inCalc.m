%% V_inCalc
function [V_in] = V_inCalc(t_vec,SIM,FLAG)

switch FLAG.InputMode
    case 1 % Step
        V_in = SIM.V_step*ones(size(t_vec));

    case 2 % Sine
        for i = 1:length(t_vec)
            t = t_vec(i);
            if t <= SIM.time_rest
                V_in(i) = SIM.V_init;
            else
                V_in(i) = SIM.V_s*sin(SIM.omega*t);
            end
        end
        
    case 3 % Ramp Step
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
    case 4 % Impulse
        for i = 1:length(t_vec)
            t = t_vec(i);
            if t <= SIM.time_rest
                V_in(i) = SIM.V_init;
            elseif t > SIM.time_rest && t <= (SIM.time_rest + SIM.time_ramp)
                V_in(i) = ((SIM.V_imp+SIM.V_init) - SIM.V_init)/SIM.time_ramp*(t-SIM.time_rest) + SIM.V_init;
            elseif t > (SIM.time_rest + SIM.time_ramp) && t <= (SIM.time_rest + SIM.Ts)
                V_in(i) = SIM.V_imp + SIM.V_init;
            elseif t > (SIM.time_rest + SIM.Ts) && t <=(SIM.time_rest + SIM.Ts + SIM.time_ramp)
                V_in(i) = (SIM.V_init - (SIM.V_imp+SIM.V_init))/SIM.time_ramp*(t-(SIM.time_rest + SIM.Ts)) + (SIM.V_imp+SIM.V_init);
            else
                V_in(i) = SIM.V_init;
            end

        end
end

end


