%% UserCurrentProfile
%
%
%


function [SIM] = UserCurrentProfile(SIM,FLAG)
%% Testing
    small = 1e-6;
    FLAG.DoXlim = 0;
    x_lim = [0,25];


%% No Noise Profile
% ---- Polarization ----
    if SIM.SimMode == 1 
        SIM.I_user_amp = SIM.C_rate*SIM.Cell_Cap*SIM.ChargeOrDischarge; % [A], Load current
        SIM.i_user_amp = SIM.I_user_amp/SIM.A_c; % [A m^-2], Load current density (flux)

        t_vec = 0;
        if SIM.initial_offset > 0 || SIM.t_ramp > 0
            amp = 0;
            if SIM.initial_offset > 0
                if SIM.t_ramp > 0
                    t_vec = [t_vec, t_vec(end)+SIM.initial_offset];
                    amp   = [amp , 0];
                else
                    t_end = t_vec(end);
                    t_vec = [t_vec, t_end + SIM.initial_offset - small, t_end + SIM.initial_offset];
                    amp   = [amp  , 0                                 , SIM.i_user_amp ];
                end
            end
            if SIM.t_ramp > 0
                t_vec = [t_vec, t_vec(end)+SIM.t_ramp];
                amp   = [amp  , SIM.i_user_amp ];
            end
        else
            amp = SIM.i_user_amp;
        end
        
        % Get Simulation Time
        if SIM.C_rate == 0
            t_final = 50; % [s], Final time
        else
            t_final = SIM.charge_frac*3600/SIM.C_rate; % [s], Sim time
        end
        SIM.profile_time_noNoise    = [t_vec, t_vec(end)+t_final];
        SIM.profile_current_noNoise = [amp  , SIM.i_user_amp    ];
    

% ---- Harmonic Perturbation ----
    elseif SIM.SimMode == 2
        SIM.I_user_amp = SIM.C_rate*SIM.Cell_Cap;
        SIM.i_user_amp = SIM.I_user_amp/SIM.A_c;

        SIM.f     = SIM.freq / (2*pi);    % [cycles s^-1],  Time to complete a period or a full sinusoid cycle
        SIM.f_s   = (2*SIM.f)*20;         % [samples s^-1], Based on Nyquist sampling theory (minimum is 2*f)
        t_sim     = SIM.N_cycles / SIM.f; % [s],            Time required for N_cycles of sinusoidal curves
        N_samples = SIM.f_s*t_sim;        % [samples],      Number of samples total
        del_time  = t_sim / N_samples;    % [s sample^-1],  Time between each sample %%%%%%% is this also 1/f_s??????????
        
        sin_time_vector = 0 : del_time : (t_sim-del_time);                     % Vector of time spaces
        sin_time_vector = sin_time_vector + SIM.initial_offset; % Add the initial offset
    
        t_vec = sin_time_vector;                        % Overall simulation vector
        amp   = SIM.i_user_amp*sin(SIM.freq*(t_vec-SIM.initial_offset)); % [A m^-2], 
        if ~SIM.initial_offset == 0
            t_vec = [0, sin_time_vector];                   % Overall simulation vector
            amp   = [0, amp];
        end
        SIM.profile_time_noNoise    = t_vec;
        SIM.profile_current_noNoise = amp;


% ---- State Space EIS ----
    elseif SIM.SimMode == 3
        SIM.I_user_amp = SIM.C_rate*SIM.Cell_Cap;
        SIM.i_user_amp = SIM.I_user_amp/SIM.A_c;

        SIM.profile_time_noNoise    = nan;
        SIM.profile_current_noNoise = nan;


% ---- Known BC Profile Controller ----
    elseif SIM.SimMode == 4
        N_steps = length(SIM.Controller_MO_File);
        for k = 1:N_steps
            switch SIM.Controller_MO_File(k).MO
                case 1 % CC
                    if SIM.Controller_MO_File(k).CorD == 'C'
                        i_user = -SIM.Controller_MO_File(k).C_rate*SIM.Cell_Cap/SIM.A_c;
                    else
                        i_user =  SIM.Controller_MO_File(k).C_rate*SIM.Cell_Cap/SIM.A_c;
                    end
                    amp    = [i_user,i_user];
                    tfinal = SIM.Controller_MO_File(k).Time_lim;
                    t_vec  = [0,tfinal];
                case 2 % CV
                    tfinal = SIM.Controller_MO_File(k).Time_lim;
                    t_vec  = [0,tfinal];
                    amp    = [nan,nan];
                case 3 % Relax
                    amp    = [0,0];
                    tfinal = SIM.Controller_MO_File(k).Time_lim;
                    t_vec  = [0,tfinal];
            end
            SIM.profile_time_noNoise_cell{k,1}    = t_vec;
            SIM.profile_current_noNoise_cell{k,1} = amp;
        end


% ---- MOO Controller ----
    elseif SIM.SimMode == 5
        SIM.profile_time_noNoise    = nan;
        SIM.profile_current_noNoise = nan;
        % Determine inside RunSimulation from MO_List

% ---- Manual Current Profile ----
    elseif SIM.SimMode == 7 % Manual Current Profile
        SIM.I_user_amp = SIM.C_rate*SIM.Cell_Cap*SIM.ChargeOrDischarge;
        SIM.i_user_amp = SIM.I_user_amp/SIM.A_c;
        SIM.i_user_OG  = SIM.i_user_amp;
        % Determine inside RunSimulation from MO_List
    

% ---- PRBS ----
    elseif SIM.SimMode == 8 
        N_PRBS   = 2^9 - 1;                  % Length of the input           [-] % New PRBS (zero-mean)
        TYPE     = 'PRBS';                   % Create a PRBS signal          [-]
        BAND     = [0 1];                    % Freq. band                    [s]
        LEVELS   = [-1 1];                   % PRBS limits                   [-]
        PRBS_gen = idinput(N_PRBS,TYPE,BAND,LEVELS);        % PRBS Demand   [A/m^2]
        C_Demand = PRBS_gen(9:SIM.PRBSLength)*SIM.PRBSAmp;  % PRBS Demand C [A/m^2] % New PRBS (zero-mean)
    
        if SIM.MakeLongPRBSSignal
            num_repeats = ceil(SIM.DesiredLength/length(C_Demand));
            C_Demand_single = C_Demand;
            for i = 1:num_repeats-1
                if mod(i,2)==1
                    C_Demand = [C_Demand ; -C_Demand_single];
                else
                    C_Demand = [C_Demand ;  C_Demand_single];
                end
            end
            C_Demand = C_Demand(1:SIM.DesiredLength);
        end
    
        t_Demand = ( ((1:length(C_Demand))-1 )*(SIM.Tswitch))' ;
        
    
        % Append a no-current entry
        if SIM.initial_offset > 0
            t_Demand = [0; t_Demand+SIM.initial_offset]; % PRBS & entry [s]
            C_Demand = [0; C_Demand];                    % PRBS & entry [A/m^2]
        end
    
        % Add Rest between some switches to help with stability
            if SIM.AddIntermediateRelaxTime
                % Find Zero-Crossings
                    change_idx = find( C_Demand(1:end-1)~=C_Demand(2:end) ); % True when the following idx doesn't match
                    change_idx = change_idx + 1; % Change occurs at the value now
    
                % Drop first change
                    change_idx = change_idx(2:end);
    
                % Find Insert index
                    mod_idx    = (SIM.NumZeroCrossingUntilNextRelax : SIM.NumZeroCrossingUntilNextRelax : length(change_idx))';
                    insert_idx = change_idx(mod_idx);
    
                % Create New Current Demand Vector
                    C_new     = C_Demand(1:insert_idx(1)-1);
                    N_changes = length(insert_idx);
                    for i = 1:N_changes-1
                        C_new = [C_new ; zeros(SIM.NumTsRelax,1) ; C_Demand(insert_idx(i):insert_idx(i+1)-1)];
                    end
                    C_new = [ C_new ; zeros(SIM.NumTsRelax,1) ; C_Demand(insert_idx(N_changes):end) ];
    
                % Create New Time Vector
                    t_new = (0:1:length(C_new)-1)*SIM.Tswitch;
            else
                C_new = C_Demand;
                t_new = t_Demand;
            end

        % Add Ramp Times
            profile_time    = t_new(1);
            profile_current = C_new(1);
        
            % Step at t^+
            for i = 2:length(t_new)
                % @ k
                profile_time(end+1,1)    = t_new(i);
                profile_current(end+1,1) = C_new(i-1);
        
                % After k
                % profile_time(end+1,1)    = t_new(i) + SIM.Tswitch * SIM.t_ramp_ratio;
                profile_time(end+1,1)    = t_new(i) + small;
                profile_current(end+1,1) = C_new(i);
            end
            profile_time    = profile_time(1:end-1);
            profile_current = profile_current(1:end-1);

        % %%%%%%%%%%%%%% Testing
        %     idx = find(profile_time<=PRBSLimit);
        %     profile_time    = profile_time(idx);
        %     profile_current = profile_current(idx);
        % %%%%%%%%%%%%%%
    
        SIM.profile_time_noNoise    = profile_time;
        SIM.profile_current_noNoise = profile_current;


% ---- EIS from Stiching PRBS ----
    elseif SIM.SimMode == 9
        SIM.profile_time_noNoise    = nan;
        SIM.profile_current_noNoise = nan;
    

% ---- EIS Ho-Kalman ----
    elseif SIM.SimMode == 10
        SIM.i_user_amp = 1;
        SIM.I_user_amp = SIM.i_user_amp * SIM.A_c;

        t_Demand    = (0:1:SIM.HK_nSamples+3) * SIM.Tsample ; % +3 is for (2 inital relax, 1 pulse)
        C_Demand    = zeros( size(t_Demand) );
        C_Demand(3) = SIM.i_user_amp;

        % Step at t^+
            profile_time    = t_Demand(1);
            profile_current = C_Demand(1);
    
            for i = 2:length(t_Demand)
                % @ k
                profile_time(end+1,1)    = t_Demand(i);
                profile_current(end+1,1) = C_Demand(i-1);
        
                % After k
                % SIM.profile_time(end+1,1)    = t_Demand(i) + SIM.Tsample * SIM.t_ramp_ratio;
                profile_time(end+1,1)    = t_Demand(i) + small;
                profile_current(end+1,1) = C_Demand(i);
            end
            profile_time    = profile_time(1:end-1);
            profile_current = profile_current(1:end-1);
    
        SIM.profile_time_noNoise    = profile_time;
        SIM.profile_current_noNoise = profile_current;
    end


%% If Adding noise (ZOH DT sim)
    if FLAG.AddInputNoise
        if SIM.SimMode == 4 % Simulations with a MO file
            [SIM.profile_time_Noisy_cell , SIM.profile_current_Noisy_cell] = getDTNoiseMO( SIM , SIM.profile_time_noNoise_cell , SIM.profile_current_noNoise_cell , small );
            SIM.profile_time_cell    = SIM.profile_time_Noisy_cell;
            SIM.profile_current_cell = SIM.profile_current_Noisy_cell;

        elseif SIM.SimMode == 8 % PRBS
            [SIM.profile_time_Noisy2D , SIM.profile_current_Noisy2D] = getDTNoisePRBS( SIM , SIM.profile_time_noNoise , SIM.profile_current_noNoise );
            SIM.profile_time_Noisy    = [];
            SIM.profile_current_Noisy = [];
            for i = 1:length(SIM.profile_time_Noisy2D)
                SIM.profile_time_Noisy    = [SIM.profile_time_Noisy   , SIM.profile_time_Noisy2D(i,:)   ];
                SIM.profile_current_Noisy = [SIM.profile_current_Noisy, SIM.profile_current_Noisy2D(i,:)];
            end
            SIM.profile_time    = SIM.profile_time_Noisy;
            SIM.profile_current = SIM.profile_current_Noisy;

        else
            [SIM.profile_time_Noisy , SIM.profile_current_Noisy] = getDTNoise( SIM , SIM.profile_time_noNoise , SIM.profile_current_noNoise , small );
            SIM.profile_time    = SIM.profile_time_Noisy;
            SIM.profile_current = SIM.profile_current_Noisy;
        end

    else
        if SIM.SimMode == 4 % Simulations with a MO file
            SIM.profile_time_cell    = SIM.profile_time_noNoise_cell;
            SIM.profile_current_cell = SIM.profile_current_noNoise_cell;

        else
            SIM.profile_time    = SIM.profile_time_noNoise;
            SIM.profile_current = SIM.profile_current_noNoise;
        end
    end


%% Plot Input Current
if FLAG.PlotUserCurrentProfiles
    if SIM.SimMode == 4 % Simulations with a MO file
        % No Noise
            t_total   = SIM.profile_time_noNoise_cell{1,1};
            amp_total = SIM.profile_current_noNoise_cell{1,1};
    
            for k = 2:N_steps
                t_end = t_total(end);
                t_next = SIM.profile_time_noNoise_cell{k,1} + t_end;
                t_total   = [t_total   , t_next];
                amp_total = [amp_total , SIM.profile_current_noNoise_cell{k,1}];
            end
            SIM.profile_time_noNoise    = t_total;
            SIM.profile_current_noNoise = amp_total;

        % Noisy
            if FLAG.AddInputNoise
                t_total   = SIM.profile_time_Noisy_cell{1,1};
                amp_total = SIM.profile_current_Noisy_cell{1,1};
        
                for k = 2:N_steps
                    t_end = t_total(end);
                    t_next = SIM.profile_time_Noisy_cell{k,1} + t_end;
                    t_total   = [t_total   ; t_next];
                    amp_total = [amp_total ; SIM.profile_current_Noisy_cell{k,1}];
                end
                SIM.profile_time_Noisy    = t_total;
                SIM.profile_current_Noisy = amp_total;
            end
    end

    % No noise
        f = figure;
        f.Position = [110   470   560   420];
        plot(SIM.profile_time_noNoise,SIM.profile_current_noNoise,'k-o','LineWidth',2,'DisplayName','No Noise')
        xlabel('Time [s]')
        ylabel('Current [A m^{-2}]')
        if FLAG.DoXlim
            xlim(x_lim)
        end
    
    % Noisy
    if FLAG.AddInputNoise
        f = figure;
        f.Position = [110+570   470   560   420];
        plot(SIM.profile_time_Noisy,SIM.profile_current_Noisy,'b--','LineWidth',2,'DisplayName','Noisy')
        xlabel('Time [s]')
        ylabel('Current [A m^{-2}]')
        if FLAG.DoXlim
            xlim(x_lim)
        end
    end
    
    % Combined
    if FLAG.AddInputNoise
        f = figure;
        f.Position = [110+2*570 470 560 420];
        hold on
        plot(SIM.profile_time_noNoise,SIM.profile_current_noNoise,'k-o','LineWidth',2,'DisplayName','No Noise')
        plot(SIM.profile_time_Noisy  ,SIM.profile_current_Noisy  ,'b--','LineWidth',2,'DisplayName','Noisy')
        xlabel('Time [s]')
        ylabel('Current [A m^{-2}]')
        lgn = legend;
        if FLAG.DoXlim
            xlim(x_lim)
        end
    end
end
end


%% getDTNoise
function [timeNoisy , currentNoisy] = getDTNoise( SIM , time_noNoise , current_noNoise , small)
    % Make DT Time vec
        t_end = time_noNoise(end);
        t_new = (0:SIM.Ts:t_end)';
        C_new   = interp1(time_noNoise , current_noNoise , t_new);
        
           % PRBS              % EIS Ho-Kalman
        if SIM.SimMode == 8 || SIM.SimMode == 10 
            t_new = t_new(1:end-1);
            C_new = C_new(2:end);
        end

    % Add noise    
        idx = find(t_new>SIM.initial_offset);
        w_k = mvnrnd(0 , SIM.Q_input , length(idx));   %w_k = w_k';
        C_new(idx) = C_new(idx) + w_k;
        % currentNoisy = amp;

    % Add t+
        timeNoisy    = t_new(1);
        currentNoisy = C_new(1);
    
        for i = 2:length(t_new)
            % @ k
            timeNoisy(end+1,1)    = t_new(i);
            currentNoisy(end+1,1) = C_new(i-1);
    
            % After k
            timeNoisy(end+1,1)    = t_new(i) + small;
            currentNoisy(end+1,1) = C_new(i);
        end
end


%% getDTNoisePRBS
function [timeNoisy , currentNoisy] = getDTNoisePRBS( SIM , time_noNoise , current_noNoise )
    % Make DT Time vec
        t_end = time_noNoise(end);
        t_new = (0:SIM.Ts:t_end)';
        C_new   = interp1(time_noNoise , current_noNoise , t_new);

        % Shift Index
            t_new = t_new(1:end-1);
            C_new = C_new(2:end);

    % Add noise    
        idx = find(t_new>=SIM.initial_offset);
        w_k = mvnrnd(0 , SIM.Q_input , length(idx));   %w_k = w_k';
        C_new(idx) = C_new(idx) + w_k;
        % currentNoisy = amp;

    % % Add t+
    %     timeNoisy    = t_new(1);
    %     currentNoisy = C_new(1);
    % 
    %     for i = 2:length(t_new)
    %         % @ k
    %         timeNoisy(end+1,1)    = t_new(i);
    %         currentNoisy(end+1,1) = C_new(i-1);
    % 
    %         % After k
    %         timeNoisy(end+1,1)    = t_new(i) + small;
    %         currentNoisy(end+1,1) = C_new(i);
    %     end

    
    % Make 2-D time and current vector
        timeNoisy    = [t_new(1:end-1) , t_new(2:end)  ];
        currentNoisy = [C_new(1:end-1) , C_new(1:end-1)];
end


%% getDTNoiseMO
function [timeNoisy , currentNoisy] = getDTNoiseMO( SIM , time_noNoise , current_noNoise , small)
    N_steps = length(SIM.Controller_MO_File);
    for k = 1:N_steps
        switch SIM.Controller_MO_File(k).MO
            case 1 % CC
                % Make DT Time vec
                    t_end = time_noNoise{k}(end);
                    t_new = (0:SIM.Ts:t_end)';
                    C_new   = interp1(time_noNoise{k} , current_noNoise{k} , t_new);

                % Add noise    
                    w_k = mvnrnd(0 , SIM.Q_input , length(t_new));   %w_k = w_k';
                    C_new = C_new + w_k;

                % Add t+
                    timeNoisy{k,1}    = t_new(1);
                    currentNoisy{k,1} = C_new(1);
                
                    for i = 2:length(t_new)
                        % @ k
                        timeNoisy{k,1}(end+1,1)    = t_new(i);
                        currentNoisy{k,1}(end+1,1) = C_new(i-1);
                
                        % After k
                        timeNoisy{k,1}(end+1,1)    = t_new(i) + small;
                        currentNoisy{k,1}(end+1,1) = C_new(i);
                    end

            case 2 % CV
                timeNoisy{k,1}    = time_noNoise{k}';
                currentNoisy{k,1} = current_noNoise{k}';
            case 3 % Relax
                timeNoisy{k,1}    = time_noNoise{k}';
                currentNoisy{k,1} = current_noNoise{k}';
        end
    end
end

%% OOOOOLLLLLLLLLDDDDDDDDDD Stuff
% Testing
    % SIM.C_rate = 20;
    % % SIM.C_rate = 0;
    % SIM.initial_offset = 10;
    % % SIM.initial_offset = 0;
    % % SIM.t_ramp = 2;
    % SIM.t_ramp = 0;
    % 
    % if SIM.C_rate == 0
    %     SIM.i_user_amp = 0;
    % else
    %     SIM.i_user_amp = 5;
    % end

    % Polarization
        % SIM.SimMode     = 1;
        % SIM.charge_frac = 1;

    % Harmonic
        % SIM.SimMode  = 2;
        % SIM.freq     = 1e-2;
        % SIM.N_cycles = 10;

    % PRBS
        % SIM.SimMode     = 8;
        % SIM.PRBSLength = 159;
        % SIM.PRBSAmp    = 1;
        % SIM.Tswitch    = 10;
        % SIM.MakeLongPRBSSignal       = 1;
        %     SIM.DesiredLength  = 4e3;
        % SIM.AddIntermediateRelaxTime = 1;
        %     SIM.NumTsRelax = 2;
        %     SIM.NumZeroCrossingUntilNextRelax = 5;
        % PRBSLimit = 200;

    % EIS Ho--Kalman
        % SIM.SimMode     = 10;
        % SIM.HK_nSamples = 800;
        % SIM.Tsample     = 1;
        % SIM.initial_offset = 0;