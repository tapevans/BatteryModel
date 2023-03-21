%% Test if PRBS is zero mean
clear all; close all; clc;

%%
power_vec = [4,5,6,7,8,9];

for PP = power_vec
%% Parameters
    Tswitch        = 1;
    PRBSLength     = 2^PP - 1;
    PRBSAmp        = 1;
    initial_offset = Tswitch;
    t_ramp_ratio   = 0.01;
    t_ramp_ratio   = Tswitch/100;


%% Make PRBS Signal
    N_PRBS  = 2^PP - 1;                 % Length of the input           [-]
    %TYPE    = 'PRBS';                   % Create a PRBS signal          [-]
    %BAND    = [0 1];                    % Freq. band                    [s]
    %LEVELS  = [-1 1];                   % PRBS limits                   [-]
    PRBS_gen= prbs(PP , N_PRBS);        % PRBS Demand   [A/m^2]

    t_Demand = ((0 : 1 : PRBSLength-1)*(Tswitch)); % PRBS Demand t [s]
    C_Demand = PRBS_gen(1 : PRBSLength)*PRBSAmp;  % PRBS Demand C [A/m^2]

    % Append a no-current entry
        if initial_offset > 0
            t_Demand = [0, t_Demand+initial_offset]; % PRBS & entry [s]
            C_Demand = [0, C_Demand];                    % PRBS & entry [A/m^2]
        end

    % Add Ramp Times
        profile_time    = t_Demand(1);
        profile_current = C_Demand(1);

    % Step at t^+
        for i = 2:length(t_Demand)
            % @ k
            profile_time(end+1,1)    = t_Demand(i);
            profile_current(end+1,1) = C_Demand(i-1);
    
            % After k
            profile_time(end+1,1)    = t_Demand(i) + Tswitch * t_ramp_ratio;
            profile_current(end+1,1) = C_Demand(i);
        end
        profile_time    = profile_time(1:end-1);
        profile_current = profile_current(1:end-1);


%% Create a discrete signal
    %t_vec  = 0 : t_ramp_ratio : profile_time(end);
    t_vec  = 0 : Tswitch/10 : profile_time(end);
    i_user = interp1(profile_time , profile_current , t_vec);


%% Integrate Signal
    coulumbs = cumtrapz( t_vec , i_user );


%% Plot Results
    figure(PP-power_vec(1)+1)
    hold on
    plot(t_vec,i_user,'k','Linewidth',2)
    yline(0,'r','Linewidth',2)
    title(['PRBS Signal ' num2str(PP) ])
    xlim([0,t_vec(end)])

    figure(PP-power_vec(1)+1+6)
    hold on
    plot(t_vec,coulumbs,'k','Linewidth',2)
    yline(0,'r','Linewidth',2)
    title(['Power ' num2str(PP) ])
    xlim([0,t_vec(end)])

end

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
