%%
% This will be used to test the Ho-Kalman algorithm
clear all; close all; clc;


%% FLAGS
% CompareMode
% 1) Impulse
% 2) Step
    FLAG.CompareMode = 1;

% HK Rank
% 1) Using Rank, r = rank(S)
% 2) All States from Hankel
% 3) Number of States
% 4) Number of short dynamics States
    FLAG.HK_rank = 4;
    
    FLAG.Use_t_minus = 0; % 1 if the impulse occurs before t_impulse; 0 if impulse occurs after
    
    FLAG.Use_Ind_ROM = 1; % 1 if using individual ROMs, 0 if using a combined ROM


%% Initial Params
    Ts          = 1; % [s], 
    TsMultiple  = 5;
    TsSim       = Ts/TsMultiple;
    N_steps_min = 100;


%% Create SS CT matrix 
    % Creation of A
    N_short_dynam = 5; 
    short_min = -2;
    short_max = -1;
    short_eig = logspace(short_min,short_max,N_short_dynam);
    
    % N_long_dynam = 10-N_short_dynam;
    % N_long_dynam = 0;
    N_long_dynam = 100;
    long_min = 1;
    long_max = 2;
    long_eig = logspace(long_min,long_max,N_long_dynam);
    
    Aeig = -1*[short_eig , long_eig];
    A = diag(Aeig);
    
    % B = zeros(length(A),1); B(1) = 1;
    B = ones(length(A),1);
    C = eye(length(A));
    %C = zeros(length(A)); C(1,1) = 1;
    D = zeros(length(C),1);
    sys_CT = ss(A,B,C,D);
    sys_DT = c2d(sys_CT,Ts);
    
    [eig_vec, eig_val] = eig(A);
    [U,S,V] = svd(A);
    
    [eig_vec_DT, eig_val_DT] = eig(sys_DT.A);
    [U_DT,S_DT,V_DT] = svd(sys_DT.A);


%% Create Time Vectors
    if 2*length(A) < N_steps_min
        t_final = N_steps_min * Ts;
    else
        t_final = 2*length(A)*Ts;
    end
    
    t_vec_ode = 0:TsSim:t_final;
    t_vec_DT  = 0:Ts:t_final;


%% Create impulse current profile
% ODE
    current = zeros(size(t_vec_ode));
    if FLAG.Use_t_minus
        idx = TsMultiple+1:1:2*TsMultiple;
    else
        idx = TsMultiple+2:1:2*TsMultiple+1;
    end
    current(idx) = ones(size(idx));
    
    if FLAG.Use_t_minus
        ImpulseCurrentProfile.t = [t_vec_ode(1:TsMultiple) , Ts-Ts/50 , t_vec_ode(TsMultiple+1:2*TsMultiple) , 2*Ts-Ts/50 , t_vec_ode(2*TsMultiple+1:end)]';
        ImpulseCurrentProfile.u = [  current(1:TsMultiple) , 0        ,   current(TsMultiple+1:2*TsMultiple) , 1          ,   current(2*TsMultiple+1:end)]';
    else
        ImpulseCurrentProfile.t = [t_vec_ode(1:TsMultiple+1) , Ts+Ts/50 , t_vec_ode(TsMultiple+2:2*TsMultiple+1) , 2*Ts+Ts/50 , t_vec_ode(2*TsMultiple+2:end)]';
        ImpulseCurrentProfile.u = [  current(1:TsMultiple+1) , 1        ,   current(TsMultiple+2:2*TsMultiple+1) , 0          ,   current(2*TsMultiple+2:end)]';
    end
    
    
    ImpulseCurrentProfile_DT.t = t_vec_DT';
    current = zeros(size(t_vec_DT)); current(2) = 1;
    ImpulseCurrentProfile_DT.u = current';


%% get Impulse response for the system
    tspan = t_vec_ode;
    SV_IC = zeros(length(A),1);
    
    % Simulation Parameters
    Tol.Abs = 1E-7;
    Tol.Rel = 1E-7;
    
    options = odeset('RelTol' ,Tol.Rel,      ...
                     'AbsTol' ,Tol.Abs,      ...
                     'MaxStep',(.99)*Ts);
    
    SOLN = ode15s(@(t,x)odeFun(t,x,sys_CT,ImpulseCurrentProfile) , tspan , SV_IC , options);
    SOLN_ODE.y = deval(SOLN,t_vec_ode);
    SOLN_ODE.x = t_vec_ode;
    SOLN_DT.y = deval(SOLN,t_vec_DT);
    SOLN_DT.x = t_vec_DT;


%% Test Plot
    if 0
        figure
        plot(SOLN.x,SOLN.y(1,:),'-ok','DisplayName','ODE')
    end


%% Test Plot
    if 0
        figure
        hold on
        plot(SOLN_ODE.x,SOLN_ODE.y(1,:),'-ok','DisplayName','ODE')
        plot(SOLN_DT.x ,SOLN_DT.y(1,:) ,'or' ,'DisplayName','DT')
        lgn = legend;
        
        figure
        hold on
        plot(SOLN_ODE.x,SOLN_ODE.y(2,:),'-ok','DisplayName','ODE')
        plot(SOLN_DT.x ,SOLN_DT.y(2,:) ,'or' ,'DisplayName','DT')
        lgn = legend;
        
        figure
        hold on
        plot(SOLN_ODE.x,SOLN_ODE.y(3,:),'-ok','DisplayName','ODE')
        plot(SOLN_DT.x ,SOLN_DT.y(3,:) ,'or' ,'DisplayName','DT')
        lgn = legend;
        
        figure
        hold on
        plot(SOLN_ODE.x,SOLN_ODE.y(end,:),'-ok','DisplayName','ODE')
        plot(SOLN_DT.x ,SOLN_DT.y(end,:) ,'or' ,'DisplayName','DT')
        lgn = legend;
    end


%% Perform Ho-Kalman
    ImpulseData_All.t = SOLN_ODE.x;
    ImpulseData_All.y = SOLN_ODE.y;
    ImpulseData = C * SOLN_DT.y;
    [sys_HK] = getHoKalmanROM(ImpulseData_All , ImpulseData , SOLN_DT.x , Ts , N_short_dynam , FLAG);
    if FLAG.Use_Ind_ROM
        for i = 1:length(sys_HK)%-1
            sys_HK_CT{i} = d2c(sys_HK{i},'tustin'); %'matched'
        end
    else
        sys_HK_CT{1} = d2c(sys_HK{end},'tustin'); %'matched'
    end


%% Analysis of HK
    if FLAG.Use_Ind_ROM
        for i = 1:length(sys_HK)%-1
            [eig_vec_HK_DT{i}, eig_val_HK_DT{i}] = eig(sys_HK{i}.A);
            [eig_vec_HK_CT{i}, eig_val_HK_CT{i}] = eig(sys_HK_CT{i}.A);
            eig_val_HK_DT_diag{i} = real(diag(eig_val_HK_DT{i}));
            eig_val_HK_CT_diag{i} = real(diag(eig_val_HK_CT{i}));
            [U_HK_DT,S_HK_DT{i},V_HK_DT] = svd(sys_HK{i}.A);
            [U_HK_CT,S_HK_CT{i},V_HK_CT] = svd(sys_HK_CT{i}.A);
        end
    else
        [eig_vec_HK_DT, eig_val_HK_DT] = eig(sys_HK{end}.A);
        [eig_vec_HK_CT, eig_val_HK_CT] = eig(sys_HK_CT{end}.A);
        eig_val_HK_DT_diag = real(diag(eig_val_HK_DT));
        eig_val_HK_CT_diag = real(diag(eig_val_HK_CT));
        [U_HK_DT,S_HK_DT,V_HK_DT] = svd(sys_HK{end}.A);
        [U_HK_CT,S_HK_CT,V_HK_CT] = svd(sys_HK_CT{end}.A);
    end


%% Get Response
    if FLAG.Use_Ind_ROM
        for i = 1:length(sys_HK)%-1
            x_red = zeros(length(sys_HK{i}.A),1);
            [z_soln_HK_DT{i} , t_soln_HK_DT{i} , x_soln_HK_DT{i}] = lsim(sys_HK{i}, ImpulseCurrentProfile_DT.u , ImpulseCurrentProfile_DT.t , x_red);
        end
    else
        x_red = zeros(length(sys_HK{end}.A),1);
        [z_soln_HK_DT{1} , t_soln_HK_DT{1} , x_soln_HK_DT{1}] = lsim(sys_HK{end}, ImpulseCurrentProfile_DT.u , ImpulseCurrentProfile_DT.t , x_red);
    end


%% Plot Response
    z_ODE = C*SOLN_ODE.y;
    for i = 1:6
        figure
        hold on
        if FLAG.Use_Ind_ROM
            plot(t_soln_HK_DT{i} , z_soln_HK_DT{i}(:,1),'or'  , 'DisplayName','HK')
        else
            plot(t_soln_HK_DT{1} , z_soln_HK_DT{1}(:,i),'or'  , 'DisplayName','HK')
        end
        plot(t_vec_ode    , z_ODE(i,:)       ,'-k' , 'DisplayName','ODE','Linewidth',2)
        xlabel('Time (s)')
        lgn = legend;
        xlim([0,50])
    end

%% Print Response
    disp('Eigenvalue Comparison')
    disp('---------------------')
    disp('Eig(A_DT)')
    disp(num2str(diag(eig_val_DT(100:105,100:105))))
    for i = 1:6
        disp(['HK System, Variable ' num2str(i)])
        disp(num2str(eig_val_HK_DT{i}))
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ODE function
function [xdot] = odeFun(t,x,sys,currentProfile)
    % get current
    u = interp1(currentProfile.t , currentProfile.u , t); 
    
    % get xdot
    xdot = sys.A * x + sys.B * u;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Ho-Kalman ROM
function [sys] = getHoKalmanROM(ImpulseData_All,  ImpulseData , t_impulse , Ts , Nrank , FLAG)  
    %% Initialize
    if FLAG.Use_t_minus
        start_idx = 1;
    else
        start_idx = 2;
    end
    [N_out,N_steps] = size(ImpulseData);


    %% Sort Data for g_k        
        z_imp    = ImpulseData;
        g_k_temp = z_imp(:,start_idx:end);
        %g_k_temp = z_imp(:,idx_vec(pulse_idx:end));    % t^-
        IC       = z_imp(:,1);
        g_k      = g_k_temp - g_k_temp(:,1);
    
        t_imp = t_impulse(start_idx:end);
        
    
    %% Plot Impulse Response
        if 0
            for i = 1:N_out
            %for i = [1,2,3,10]
                figure
                hold on
                plot(ImpulseData_All.t , ImpulseData_All.y(i,:) ,'-ok'                      ,'DisplayName','ODE')
                plot(t_imp             , g_k(i,:)+IC(i,1)       ,'or' ,'MarkerFaceColor','r','DisplayName','Selected Data')
                %plot(t_imp(idx_vec(pulse_idx:end),:)  ,g_k(i,:)+IC(i,1),'or','MarkerFaceColor','r') % t^-
                title(['Data Used in Ho-Kalman Reduction, State ' num2str(i)])
                lgn = legend;
                lgn.Location = 'best';
            end
        end

if FLAG.Use_Ind_ROM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get ROM for each desired outputs separately (Combined with Voltage)
    if N_out > Nrank
        loop_end = Nrank+1;
    else
        loop_end = Nrank;
    end

        for OO = 1:loop_end
            %% Get Hankle Matrix
            H = myHankle(g_k(OO,:));
    
    
            %% Convert to SS DT sys
            % SVD of Hankle
            [Nrows, Ncolms] = size(H);
            N_outputs = Nrows/Ncolms;
            N_in = 1; %%% !!! Hardcoded but always true for batteries
    
            [U,S,V] = svd(H);
            switch FLAG.HK_rank
                case 1
                    r = rank(S);
                case 2
                    r = Ncolms;
                case 3
                    r = N_out;
                case 4
                    r = Nrank;
            end
            r = 2;
            r

            %singVals(OO) = S(r,r);
    
            U_colm = U(:,1:r);
            S_colm = S(1:r,1:r);
            V_colm = V(:,1:r);
    
            % Balanced Reduction
            S_sqrt = S_colm.^0.5;
            obsv_r = U_colm*S_sqrt;
            cont_r = S_sqrt*V_colm';
    
            % C_r
            C_r = obsv_r(1:N_outputs,:   );
    
            % B_r
            B_r = cont_r(:        ,N_in);
    
            % A_r
            PP   = obsv_r(1:end-N_outputs,:);
            PP_p = obsv_r(N_outputs+1:end,:);
            threshold = 1e-7;
            PP_inv = pinv(PP,threshold);
            A_r = PP_inv*PP_p;
    
            % D_r ~ !!!!! I'm assuming this is correct
            D_r = zeros(N_outputs,N_in);

%             if FLAG.PlotSingVal
%                 figure
%                 semilogy(1:1:length(diag(S)) , diag(S),'-k','Linewidth',2)
%                 xlim([0,length(diag(S))])
%                 ylabel('Singular Values')
%                 title([RESULTS.Labels.title{OO} ' Singular Values of \Sigma for Hankel'])
%             end
    
    
            %% Return the system
            sys{OO} = ss(A_r, B_r, C_r, D_r, Ts);
        end

else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Get ROM for all desired outputs combined
            %% Get Hankle Matrix
            H = myHankle(g_k);
    
    
            %% Convert to SS DT sys
            % SVD of Hankle
            [Nrows, Ncolms] = size(H);
            N_outputs = Nrows/Ncolms;
            N_in = 1; %%% !!! Hardcoded but always true for batteries
    
            [U,S,V] = svd(H);
            switch FLAG.HK_rank
                case 1
                    r = rank(S);
                case 2
                    r = Ncolms;
                case 3
                    r = N_out;
                case 4
                    r = Nrank;
            end
            r
            
            %singVals = S(r,r)
    
            U_colm = U(:,1:r);
            S_colm = S(1:r,1:r);
            V_colm = V(:,1:r);
    
            % Balanced Reduction
            S_sqrt = S_colm.^0.5;
            obsv_r = U_colm*S_sqrt;
            cont_r = S_sqrt*V_colm';
    
            % C_r
            C_r = obsv_r(1:N_outputs,:   );
    
            % B_r
            B_r = cont_r(:        ,N_in);
    
            % A_r
            PP   = obsv_r(1:end-N_outputs,:);
            PP_p = obsv_r(N_outputs+1:end,:);
            threshold = 1e-7;
            PP_inv = pinv(PP,threshold);
            A_r = PP_inv*PP_p;
    
            % D_r ~ !!!!! I'm assuming this is correct
            D_r = zeros(N_outputs,N_in);
    
            if 0
                figure
                semilogy(1:1:length(diag(S)) , diag(S),'-k','Linewidth',2)
                xlim([0,length(diag(S))])
                ylabel('Singular Values')
                title('All ROM Singular Values of \Sigma for Hankel')
            end
    
    
            %% Return the system
            sys{1} = ss(A_r, B_r, C_r, D_r, Ts);
%             sys{1 , end+1} = ss(A_r, B_r, C_r, D_r, Ts);
%             sys = ss(A_r, B_r, C_r, D_r, Ts);

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate SS
%% Create step current profile
%% get Step response