function [u , z] = ROMplantData(SIM,FLAG,N,P,RESULTS)
%% Get Noise
    [InputSig_w , InputSig_v] = getSlinkNoise(SIM,N,P,FLAG);

    % Find times to obtain noise
    idx = [1 4:3:length(InputSig_w)];
    w_k = InputSig_w(idx,2)';
    v_k = InputSig_v(idx,2)';

    %w_k = reshape(w_k,1,[]);
    %v_k = reshape(v_k,1,[]);


%% Get u_k
    u_k = SIM.InputSignal(:,2);
    u_k = reshape(u_k,1,[]);


%% Get ROM
    [sys_HK,~] = getHoKalmanROM(SIM,N,P,FLAG,RESULTS);
    A_DT = sys_HK{1}.A;
    B_DT = sys_HK{1}.B;
    C_DT = sys_HK{1}.C(1,:);


%% Initialize
    N_steps = length(u_k);
    [N_measur, N_states] = size(C_DT);
    
    x_EST_sys = nan(N_states,N_steps);
    z_k       = nan(N_measur,N_steps);


%% Actuate System
    for i = 2:N_steps
        x_EST_sys(:,i) = A_DT * x_EST_sys(:,i-1)  +  B_DT * u_k(i-1)  +  w_k(:,i-1);
        z_k(:,i)       = C_DT * x_EST_sys(:,i)  +  v_k(:,i);
    end


%% Save Outputs
    u = u_k;
    z = z_k;

    
end