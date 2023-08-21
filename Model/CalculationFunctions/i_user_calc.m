function i_user = i_user_calc(t_in,SIM,FLAG,i_user_input)
% Use this function to update i_user inside the ode function

%    % Polarization      % PRBS              % EIS from Stitching PRBS     % EIS Ho-Kalman
% if SIM.SimMode == 1 || SIM.SimMode == 8 || SIM.SimMode == 9           || SIM.SimMode == 10
   % Polarization      % EIS from Stitching PRBS     % EIS Ho-Kalman
if SIM.SimMode == 1 || SIM.SimMode == 9           || SIM.SimMode == 10
    i_user = interp1(SIM.profile_time , SIM.profile_current , t_in);


elseif SIM.SimMode == 2 % ---- Harmonic Perturbation ----
    if FLAG.AddInputNoise
        i_user = interp1(SIM.profile_time , SIM.profile_current , t_in);
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


elseif SIM.SimMode == 4 % ---- Known BC Profile Controller ----
    if SIM.Controller_MO_File(SIM.current_MO_step).MO == 2 % CV
        i_user = i_user_input; %%%%%% How will this work for vector of input values?
    else
        i_user = interp1(SIM.profile_time , SIM.profile_current , t_in);
    end


elseif SIM.SimMode == 7 % ---- Manual Current Profile ----
    if isfield(SIM,'profile_time')
        i_user = interp1(SIM.profile_time , SIM.profile_current , t_in);
    else
        i_user = 0; %%%%%%%%% Don't know if this will work
    end


elseif SIM.SimMode == 8 % ---- PRBS ----
    i_user = i_user_input; %%%%%% How will this work for vector of input values?


else
    %%%%%%%% Don't know when this condition occurs
    i_user = SIM.Amp * ones(1,length(t_in));
    disp('Hit the else condition in i_user')
end

end

%% OOOOOOOLLLLLLLLLDDDDDDDDDD

% % i_user = zeros(1,length(t_in));
% for j = 1:length(t_in)
%     % t = t_in(j);
% 
%     % ---- Polarization ----
%     if SIM.SimMode == 1 
%     %     if t < SIM.initial_offset % Zero initial then polarization
%     %         i_user(j) = 0;
%     %     elseif ( t>=SIM.initial_offset ) && ( t< (SIM.initial_offset+SIM.t_ramp) )
%     %         i_user(j) = (SIM.i_user_amp / SIM.t_ramp) * t;
%     %     else
%     %         i_user(j) = SIM.i_user_amp;
%     %     end    
% 
%     % ---- Harmonic Perturbation ----
%     elseif SIM.SimMode == 2 % Zero initial then sinusoidal (EIS)
%         % if t < SIM.initial_offset
%         %     i_user(j) = 0;
%         % else
%         %     i_user(j) = SIM.i_user_amp*sin(SIM.freq*(t-SIM.initial_offset)); % [A m^-2], 
%         % end
% 
%     % ---- Manual Current Profile ----
%     elseif SIM.SimMode == 7 
%         % if isfield(SIM,'profile_time')
%         %     i_user(j) = interp1(SIM.profile_time , SIM.profile_current , t);
%         % else
%         %     i_user(j) = 0;
%         % end
% 
%     % ---- PRBS ----
%     elseif SIM.SimMode == 8 
%         % i_user(j) = interp1(SIM.profile_time , SIM.profile_current , t);
% 
%     % ---- EIS from Stitching PRBS ---- 
%     elseif SIM.SimMode == 9 
%         % i_user(j) = interp1(SIM.profile_time , SIM.profile_current , t);
% 
%     % ---- EIS Ho-Kalman ----
%     elseif SIM.SimMode == 10 
%         % i_user(j) = interp1(SIM.profile_time , SIM.profile_current , t);
% 
%     else
%         i_user(j) = SIM.Amp;
%     end
% end