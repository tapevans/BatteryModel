function [t, u , x , z_all] = ROMplantData(SIM,FLAG,N,P,RESULTS)
    %% Get Noise
        [InputSig_w , InputSig_v] = getSlinkNoise(SIM,N,P,FLAG);
    
    % Find times to obtain noise
        idx = [1 4:3:length(InputSig_w)];
        w_k = InputSig_w(idx,2)';
        v_k = InputSig_v(idx,2)';

%         idx = [1 3:3:length(InputSig_w)];
%         t   = InputSig_w(idx,1)';
    
    
    %% Get u_k
        t   = SIM.InputSignal(:,1);
        u_k = SIM.InputSignal(:,2);
        u_k = reshape(u_k,1,[]);
        
        u = u_k;
        N_steps = length(u_k);

        z_init = SIM.y_0_FOM;
    
    
    %% Get ROM
    [sys_HK,~] = getHoKalmanROM(SIM,N,P,FLAG,RESULTS);
    if FLAG.EST.SepHK
        for j = 1:N.DesOut
            %% Get SS Matricies
            A_DT = sys_HK{j}.A;
            B_DT = sys_HK{j}.B;
            C_DT = sys_HK{j}.C(1,:);
    
    
            %% Initialize
            [N_measur, N_states] = size(C_DT);
    
            x_EST_sys = nan(N_states,N_steps);
            z_k       = nan(N_measur,N_steps);

            x_EST_sys(:,1) = zeros(N_states,1);
            z_k(:,1)       = z_init(j,1);

    
            %% Actuate System
            for i = 2:N_steps
                x_EST_sys(:,i) = A_DT * x_EST_sys(:,i-1)  +  B_DT * u_k(i-1)  +  w_k(:,i-1);
                z_k(:,i)       = C_DT * x_EST_sys(:,i)                        +  v_k(:,i);
            end
    
    
            %% Save Outputs
            x{j}     = x_EST_sys;
            z_all{j} = sys_HK{j}.C * x_EST_sys + z_init([1,j],1);
    
        end
    
    else
        %% Get SS Matricies
        A_DT = sys_HK{end}.A;
        B_DT = sys_HK{end}.B;
        C_DT = sys_HK{end}.C(1,:);
    
    
        %% Initialize
        N_steps = length(u_k);
        [N_measur, N_states] = size(C_DT);
    
        x_EST_sys = nan(N_states,N_steps);
        z_k       = nan(N_measur,N_steps);
    
    
        %% Actuate System
        for i = 2:N_steps
            x_EST_sys(:,i) = A_DT * x_EST_sys(:,i-1)  +  B_DT * u_k(i-1)  +  w_k(:,i-1);
            z_k(:,i)       = C_DT * x_EST_sys(:,i)                        +  v_k(:,i);
        end
    
    
        %% Save Outputs
        x{1}     = x_EST_sys;
        z_all{1} = sys_HK{end}.C * x_EST_sys;
        
    
    end
end