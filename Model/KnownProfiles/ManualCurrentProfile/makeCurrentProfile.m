%% Make Current Profile
% Description: 
% 
% Inputs:
%   *
%   *
%   *
%
% Outputs:
%   *
%   *
        
function [profile_save_filepath] = makeCurrentProfile(profile_save_filepath, SIM)
    %% Check if this profile has already been created
        %% Fix file name
        if SIM.plating_refine
            new_name = profile_save_filepath(1:end-4); % removes '.mat'
            new_name = [new_name '_CurrentProfile_Input.mat'];
            profile_save_filepath = new_name;
        else
            % Nothing here
        end
        %% Check if it exists
        profile_exists = isfile(profile_save_filepath);
        if profile_exists
            %%%%% Later add something that allows me to overwrite existing file            
        else
            if SIM.plating_refine
                % Region is used more as the piecewise function
                
                % Profile includes the intermediate time steps for ramping
                % since ode can't handle sharp changes i_user. This is used
                % for interpolation
                
                t_final = 3600;
                % Region Time 
                    region_time_vec = linspace(0,t_final, SIM.N_regions+1);
                % Region Current
                    region_current_vec = SIM.i_user_amp * ones(SIM.N_regions,1);
                % Profile Time
                    profile_time = 0;
                    i = 1;
                    profile_time(end+1) = region_time_vec(i) + SIM.ramp_time/2;

                    for i = 2:SIM.N_regions
                        profile_time(end+1) = region_time_vec(i) - SIM.ramp_time/2;
                        profile_time(end+1) = region_time_vec(i);
                        profile_time(end+1) = region_time_vec(i) + SIM.ramp_time/2;
                    end

                    i = SIM.N_regions+1;
                    profile_time(end+1) = region_time_vec(i);
                % Profile Current
                    profile_current = SIM.i_user_amp * ones(size(profile_time));
                    profile_current(1) = 0;
            else
                % Whatever I want to do here to make a new current profile
            end
            clear SIM
            SIM.SimMode = 0;
            save(profile_save_filepath,'region_time_vec','region_current_vec','profile_time','profile_current','t_final','profile_save_filepath','SIM')
        end

end

% SIM.region_time_vec        = region_time_vec;
% SIM.region_current_vec     = region_current_vec;
% SIM.profile_time           = profile_time;
% SIM.profile_current        = profile_current;
% SIM.tspan                  = [0, t_final];
% SIM.Input_Profile_filepath = profile_save_filepath;
