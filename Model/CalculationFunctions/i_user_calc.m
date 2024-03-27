function i_user = i_user_calc(t_in, i_user_input, SIM__SimMode, SIM__profile_time, SIM__profile_current, SIM__Amp, SIM, FLAG)
% Use this function to update i_user inside the ode function

%    % Polarization      % PRBS              % EIS from Stitching PRBS     % EIS Ho-Kalman
% if SIM.SimMode == 1 || SIM.SimMode == 8 || SIM.SimMode == 9           || SIM.SimMode == 10
   % Polarization      % EIS from Stitching PRBS     % EIS Ho-Kalman       % MultiLevelInput
if SIM__SimMode == 1 || SIM__SimMode == 9           || SIM__SimMode == 10 || SIM__SimMode == 11
    i_user = interp1(SIM__profile_time , SIM__profile_current , t_in);


elseif SIM__SimMode == 2 % ---- Harmonic Perturbation ----
    if FLAG.AddInputNoise
        i_user = interp1(SIM__profile_time , SIM__profile_current , t_in);
    else
        i_user = zeros(1,length(t_in));
        for j = 1:length(t_in)
            t = t_in(j);
            if t < SIM.initial_offset
                i_user(j) = 0;
            else
                i_user(j) = SIM.i_user_amp*sin(SIM.freq*(t-SIM.initial_offset)); % [A m^-2], 
            end
        end
    end


elseif SIM__SimMode == 4 % ---- Known BC Profile Controller ----
    if SIM.Controller_MO_File(SIM.current_MO_step).MO == 2 % CV
        i_user = i_user_input; %%%%%% How will this work for vector of input values?
    else
        i_user = interp1(SIM__profile_time , SIM__profile_current , t_in);
    end


elseif SIM__SimMode == 7 % ---- Manual Current Profile ----
    if isfield(SIM,'profile_time')
        i_user = interp1(SIM__profile_time , SIM__profile_current , t_in);
    else
        i_user = 0; %%%%%%%%% Don't know if this will work
    end


elseif SIM__SimMode == 8 % ---- PRBS ----
    i_user = i_user_input; %%%%%% How will this work for vector of input values?


else
    %%%%%%%% Don't know when this condition occurs
    i_user = SIM__Amp * ones(1,length(t_in));
    disp('Hit the else condition in i_user')
end

end