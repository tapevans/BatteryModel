function i_user = i_user_calc(t_in,SIM)
% Use this function to update i_user inside the ode function
i_user = zeros(1,length(t_in));
for j = 1:length(t_in)
    t = t_in(j);

    % ---- Polarization ----
    if SIM.SimMode == 1 
        if t < SIM.initial_offset % Zero initial then polarization
            i_user(j) = 0;
        elseif ( t>=SIM.initial_offset ) && ( t< (SIM.initial_offset+SIM.t_ramp) )
            i_user(j) = (SIM.i_user_amp / SIM.t_ramp) * t;
        else
            i_user(j) = SIM.i_user_amp;
        end    

    % ---- Harmonic Perturbation ----
    elseif SIM.SimMode == 2 % Zero initial then sinusoidal (EIS)
        if t < SIM.initial_offset
            i_user(j) = 0;
        else
            i_user(j) = SIM.i_user_amp*sin(SIM.freq*(t-SIM.initial_offset)); % [A m^-2], 
        end

    % ---- Manual Current Profile ----
    elseif SIM.SimMode == 7 
        if isfield(SIM,'profile_time')
            i_user(j) = interp1(SIM.profile_time , SIM.profile_current , t);
        else
            i_user(j) = 0;
        end

    % ---- PRBS ----
    elseif SIM.SimMode == 8 
        i_user(j) = interp1(SIM.profile_time , SIM.profile_current , t);

    % ---- EIS from Stiching PRBS ---- 
    elseif SIM.SimMode == 9 
        i_user(j) = interp1(SIM.profile_time , SIM.profile_current , t);

    % ---- EIS Ho-Kalman ----
    elseif SIM.SimMode == 10 
        i_user(j) = interp1(SIM.profile_time , SIM.profile_current , t);

    else
        i_user(j) = SIM.Amp;
    end
end

end