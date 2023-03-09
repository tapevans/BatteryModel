%% Get Ho-Kalman ROM
function [sys,singVals] = getHoKalmanROM(SIM,N,P,FLAG,RESULTS)
%% Load Impulse Data
    filename = getImpulseFilename(FLAG);
    simsys = load(filename);


%% Sort Data for g_k
    % Extract g_k
    idx_vec = 1:simsys.SIM.TsMultiple:length(simsys.t_soln);

    i_user_DT = simsys.i_user(idx_vec);
    pulse_idx = find(i_user_DT ~= 0 );
    
    z_imp    = SIM.OutputMatrix * simsys.SV_soln';
    g_k_temp = z_imp(:,idx_vec(pulse_idx-1:end));
    %g_k_temp = z_imp(:,idx_vec(pulse_idx:end));    % t^-
    IC       = z_imp(:,1);
    g_k      = g_k_temp - g_k_temp(:,1);

    t_imp = simsys.t_soln;


%% Plot Impulse Response
    if FLAG.Analysis.PlotImp
        for i = 1:N.DesOut
        figure
        hold on
        plot(t_imp                            ,z_imp(i,:)      ,'-ok','DisplayName','ODE')
        plot(t_imp(idx_vec(pulse_idx-1:end),:),g_k(i,:)+IC(i,1),'or','MarkerFaceColor','r','DisplayName','Selected Data')
        %plot(t_imp(idx_vec(pulse_idx:end),:)  ,g_k(i,:)+IC(i,1),'or','MarkerFaceColor','r') % t^-
        title(['Data Used in Ho-Kalman Reduction ' RESULTS.Labels.title{i}])
        lgn = legend;
        lgn.Location = 'best';
        end
    end


%% Initialize Function Outputs
    sys      = cell(1,N.DesOut+1);
    singVals = nan(1, N.DesOut+1);

%% Optimal
if FLAG.UseOptimal
    [ROMError_filename] = getROMErrorFilename(FLAG);
    if isfile(ROMError_filename)
        optimal_r_obj = load(ROMError_filename,'optimal_r');
    else
        optimal_r_obj = 0;
    end
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get ROM for each desired outputs separately (Combined with Voltage)
        for OO = 1:N.DesOut
            %% Get Hankle Matrix
            H = myHankle(g_k([1,OO],:));
    
    
            %% Convert to SS DT sys
            % SVD of Hankle
            [Nrows, Ncolms] = size(H);
            N_outputs = Nrows/Ncolms;
            N_in = 1; %%% !!! Hardcoded but always true for batteries
    
            [U,S,V] = svd(H);
            r = rank(S);
            %r = rank(S,1.15e-7);
%             if FLAG.UseInput_r
%                 r = SIM.Input_r;
%             else
%                 if FLAG.UseOptimal
%                     r = optimal_r_obj.optimal_r.ind.T100(OO);
%                 else
%                     %r = 38;
%                     r = 18;
%                 end
%             end
%             if OO == P.delta_C_Li
%                 r = 18;
%             end
%             % Delete Later!!!!!!!!!!!!!!!!!!!!
%             switch OO
%                 case 1
%                     r = 23;
%                 case 2
%                     r = 18;
%                 case 3
%                     r = 29;
%                 case 4
%                     r = 14;
%                 case 5
%                     r = 49;
%                 case 6
%                     r = 25;
%             end



            singVals(OO) = S(r,r);
    
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

            if FLAG.PlotSingVal
                figure
                semilogy(1:1:length(diag(S)) , diag(S),'-k','Linewidth',2)
                xlim([0,length(diag(S))])
                ylabel('Singular Values')
                title([RESULTS.Labels.title{OO} ' Singular Values of \Sigma for Hankel'])
            end
    
    
            %% Return the system
            sys{OO} = ss(A_r, B_r, C_r, D_r, SIM.Ts);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get ROM for all desired outputs combined
        %% Get Hankle Matrix
        H = myHankle(g_k);


        %% Convert to SS DT sys
        % SVD of Hankle
        [Nrows, Ncolms] = size(H);
        N_outputs = Nrows/Ncolms;
        N_in = 1; %%% !!! Hardcoded but always true for batteries

        [U,S,V] = svd(H);
        r = rank(S);
        %r = rank(S,1.15e-7);
%         if FLAG.UseInput_r
%             r = SIM.Input_r;
%         else
%             if FLAG.UseOptimal
%                 r = optimal_r_obj.optimal_r.ind.AllSims(end);
%             else
%                 %r = 38;
%                 %r = 18;
%                 r = 9;
%             end
%         end
        
        singVals(end) = S(r,r);

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

        if FLAG.PlotSingVal
            figure
            semilogy(1:1:length(diag(S)) , diag(S),'-k','Linewidth',2)
            xlim([0,length(diag(S))])
            ylabel('Singular Values')
            title('All ROM Singular Values of \Sigma for Hankel')
        end


        %% Return the system
        sys{1 , N.DesOut+1} = ss(A_r, B_r, C_r, D_r, SIM.Ts);

end

