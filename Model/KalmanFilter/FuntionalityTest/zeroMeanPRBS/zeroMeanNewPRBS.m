%% Test if PRBS is zero mean
clear all; close all; clc;

%% Parameters
    Tswitch        = 1;
    %PRBSLength     = 159;
    PRBSLength     = 77;
    PRBSAmp        = 1;
    initial_offset = Tswitch;
    t_ramp_ratio   = 1/50;

    AddIntermediateRelaxTime = 1;
        NumTsRelax = 2;
        NumZeroCrossingUntilNextRelax = 5;


%% Make PRBS Signal
    N_PRBS  = 2^9 - 1;                 % Length of the input           [-]
    TYPE    = 'PRBS';                   % Create a PRBS signal          [-]
    BAND    = [0 1];                    % Freq. band                    [s]
    LEVELS  = [-1 1];                   % PRBS limits                   [-]
    PRBS_gen= idinput(N_PRBS,TYPE,BAND,LEVELS);         % PRBS Demand   [A/m^2]

    %C_Demand = PRBS_gen(1 : PRBSLength)*PRBSAmp;  % PRBS Demand C [A/m^2]
    %t_Demand = ((0 : 1 : PRBSLength-1)*(Tswitch))'; % PRBS Demand t [s]
    C_Demand = PRBS_gen(9 : PRBSLength)*PRBSAmp;  % PRBS Demand C [A/m^2]
    t_Demand = ((0 : 1 : PRBSLength-9)*(Tswitch))'; % PRBS Demand t [s]
    

    %% Append a no-current entry
        if initial_offset > 0
            t_Demand = [0; t_Demand+initial_offset]; % PRBS & entry [s]
            C_Demand = [0; C_Demand];                % PRBS & entry [A/m^2]
        end
    
    %% Add Rest
        if AddIntermediateRelaxTime
            % Find Zero-Crossings
            change_idx = find( C_Demand(1:end-1)~=C_Demand(2:end) ); % True when the following idx doesn't match
            change_idx = change_idx + 1; % Change occurs at the value now
            % Drop first change
            change_idx = change_idx(2:end);
            % Find Insert index
            mod_idx = (NumZeroCrossingUntilNextRelax : NumZeroCrossingUntilNextRelax : length(change_idx))';
            insert_idx = change_idx(mod_idx);
            % Create New Current Demand Vector
            C_new = C_Demand(1:insert_idx(1)-1);
            N_changes = length(insert_idx);
            for i = 1:N_changes-1
                C_new = [C_new ; zeros(NumTsRelax,1) ; C_Demand(insert_idx(i):insert_idx(i+1)-1)];
            end
            C_new = [ C_new ; zeros(NumTsRelax,1) ; C_Demand(insert_idx(N_changes):end) ];
            % Create New Time Vector
            t_new = (0:1:length(C_new)-1)*Tswitch;
        else
            C_new = C_Demand;
            t_new = t_Demand;
        end


    %% Add Ramp Times
        profile_time    = t_new(1);
        profile_current = C_new(1);

    %% Step at t^+
        for i = 2:length(t_new)
            % @ k
            profile_time(end+1,1)    = t_new(i);
            profile_current(end+1,1) = C_new(i-1);
    
            % After k
            profile_time(end+1,1)    = t_new(i) + Tswitch * t_ramp_ratio;
            profile_current(end+1,1) = C_new(i);
        end
        profile_time    = profile_time(1:end-1);
        profile_current = profile_current(1:end-1);


%% Create a discrete signal
    %t_vec  = 0 : t_ramp_ratio : profile_time(end);
    t_vec  = 0 : Tswitch/10 : profile_time(end);
    i_user = interp1(profile_time , profile_current , t_vec);


%% Integrate Signal
    coulombs = cumtrapz( t_vec , i_user );


%% Plot Results
    figure
    hold on
    plot(t_vec,i_user,'k','Linewidth',2)
    %yline(0,'r','Linewidth',2)
    title(['PRBS Signal'])
    xlim([0,t_vec(end)])

    figure
    hold on
    plot(t_vec,coulombs,'k','Linewidth',2)
    yline(0,'r','Linewidth',2)
    title(['Cumulative Coulombs'])
    xlim([0,t_vec(end)])


%% Arrange Figures
FigArrange = 1;
if FigArrange == 1
    fig = gcf;
    NumFig = fig.Number;
    if NumFig == 1 && isempty(fig.Children)
        close(fig)
    else
        Ncol = 3;
        
        for i = 1:NumFig
            f = figure(i);
            k = mod(i-1,Ncol);
            row = mod(fix((i-1)/Ncol),2);
            if row == 0
                r = 575;
    %             r = 540;
            elseif row == 1
                r = 62;
            end
            f.Position = [k*575+15 r 560 420];
        end
    end
end
