%% CC_CV_CC_Controller
% Modes of Operation (MO)
% 0) Simulation Ended
% 1) Constant Current
% 2) Constant Voltage
% 3) Relaxation
function [i_user] = Controller_CC_CV_CC(t, SV, AN , SEP, CA , EL , P , N , SIM , CONS , FLAG , PROPS)
% Parameters
	K_p = 0.5;
	K_i = 0.5;
	K_d = 0.5;

%% Setup Variables used in tracking 
	persistent current_MO_step
	if isempty(current_MO_step)
		current_MO_step = 1;
	end

	persistent time_history
	if isempty(time_history)
		time_history = zeros(10,1);
	end

	persistent error_history
	if isempty(error_history)
		error_history = zeros(10,1);
	end

	persistent i_user_history
	if isempty(i_user_history)
		i_user_history = zeros(10,1);
	end

	persistent MO_StartEnd_Time
	if isempty(MO_StartEnd_Time)
		MO_StartEnd_Time = zeros(length(SIM.Controller_MO_File),2);
	end

%% Cell Voltage Calc
	cell_voltage = SV(P.phi_ed,end) - SV(P.phi_ed,1);

%% Check if end conditions for current MO have been reached
	MO = SIM.Controller_MO_File(current_MO_step).MO;
	end_reached = 0;
	if MO == 1 % CC
		% Check Voltage Limit
			if SIM.Controller_MO_File(current_MO_step).CorD == 'C'
				if cell_voltage > SIM.Controller_MO_File(current_MO_step).Volt_lim
					end_reached = 1;
				end
			elseif SIM.Controller_MO_File(current_MO_step).CorD == 'D'
				if cell_voltage < SIM.Controller_MO_File(current_MO_step).Volt_lim
					end_reached = 1;
				end	
			end
		% Check Time Limit
			start_time = MO_StartEnd_Time(current_MO_step,1);
			delta_time = t - start_time;
			if delta_time > SIM.Controller_MO_File(current_MO_step).Time_lim
				end_reached = 1;
			end	
	elseif MO == 2 % CV
		% Check Time Limit
			start_time = MO_StartEnd_Time(current_MO_step,1);
			delta_time = t - start_time;
			if delta_time > SIM.Controller_MO_File(current_MO_step).Time_lim
				end_reached = 1;
			end	
		% Check Change in Variables Tolerance
			% !!!!!Not Implemented
	elseif MO == 3 % Relaxation 
		% No end condition at this time. It is used to run out the last of
		% tspan for the ode
	else
		disp('do not know what mode of operation I am in')
	end

% If the end condition has been reached, move to the next MO
if end_reached
    MO_StartEnd_Time(current_MO_step,2) = t;
    if current_MO_step ~= length(SIM.Controller_MO_File)
        current_MO_step = current_MO_step + 1;
        MO_StartEnd_Time(current_MO_step,1) = t;
        MO = SIM.Controller_MO_File(current_MO_step).MO;
    else
        % !!!!!!!!Figure out how to break out of the ODE
            % Currently using Relaxation to finish sim
    end
end


%% Error Calculation
    error = 0;
	if MO == 1 % CC
		% There is no error that needs to be calculated, just out the proper current
		error = 0;
	elseif MO == 2 % CV
		voltage_ref = SIM.Controller_MO_File(current_MO_step).Volt_ref;
		error = voltage_ref - cell_voltage;
	elseif MO == 3 % Relaxation
		error = 0;
	else
		disp('do not know what mode of operation I am in again')
	end

	error_history  = [error_history(2:end) ; error];
	time_history   = [time_history(2:end)  ; t      ];


%% Output Current Calculation
if MO == 2
    MO;
end
	if MO == 1 % CC
		% There is no error that needs to be calculated, just out the proper current
        if SIM.Controller_MO_File(current_MO_step).CorD == 'C'
            i_user = -SIM.Controller_MO_File(current_MO_step).C_rate*SIM.Cell_Cap/SIM.A_c;
        else
            i_user =  SIM.Controller_MO_File(current_MO_step).C_rate*SIM.Cell_Cap/SIM.A_c;
        end
	elseif MO == 2 % CV
        new_del_t = abs(time_history(end) - time_history(end -1));
        if new_del_t <= 1e-16 % handle dividing by zero in the denominator
            percent_change_i_user =   K_p * error ...
                                    + K_i * trapz(time_history,error_history);            
        else
            percent_change_i_user =   K_p * error ...
                                    + K_i * trapz(time_history,error_history) ...
                                    + K_d * (error_history(end) - error_history(end -1))/(time_history(end) - time_history(end -1));            
        end
%         i_user = i_user_history(end) + percent_change_i_user * i_user_history(end);
        i_user = i_user_history(end) - percent_change_i_user;
	elseif MO == 3 % Relaxation
		i_user = 0;
	else
		disp('do not know what mode of operation I am in AGAIN')
	end
	
	i_user_history =  [i_user_history(2:end) ; i_user];
	
%% Troubleshooting
% Write stuff to a file
    fileID = fopen(SIM.TroubleshootFilename,'a+');
    fprintf(fileID,'%12e  %6i  %6i  %8.4f  %8.4f  %12e\n',t     ,current_MO_step, MO    , cell_voltage , i_user   , error  );
    fclose(fileID);


end