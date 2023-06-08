%% Get Ho-Kalman ROM
function [sys,singVals] = getHoKalman(t_soln , i_user , z_imp , r , SIM , FLAG)
%% Sort Data for g_k
% Extract g_k
    idx_vec = 1 : SIM.TsMultiple : length(t_soln);
    
    i_user_DT = i_user(idx_vec);
    pulse_idx = find(i_user_DT ~= 0 );
    
    g_k_temp = z_imp(:,idx_vec(pulse_idx-1:end));
    IC       = z_imp(:,1);
    g_k      = g_k_temp - g_k_temp(:,1);


%% Adjust for Order of Magnitude
    range    = (max(g_k')-min(g_k'));
    g_k      = inv(diag(range))*g_k;
    multiple = diag(range);


%% Get Hankel Matrix
    H = myHankel(g_k);


%% Convert to SS DT sys
% SVD of Hankel
    [Nrows, Ncolms] = size(H);
    N_outputs = Nrows/Ncolms;
    N_in = 1; %%% !!! Hardcoded but always true for batteries
    
    [U,S,V] = svd(H);
    
    if r > Ncolms % Asking for a higher order model than possible from the data
        r = Ncolms;
    end

    singVals = S(r,r);

    U_colm = U(:,1:r);
    S_colm = S(1:r,1:r);
    V_colm = V(:,1:r);
    
% Balanced Reduction
    S_sqrt = S_colm.^0.5;
    obsv_r = U_colm*S_sqrt;
    cont_r = S_sqrt*V_colm';


%% Get State-Space matrices
% C_r
    C_r = obsv_r( 1 , : );
    C_r = multiple * C_r;
    
% B_r
    B_r = cont_r(: ,N_in);
    
% A_r
    PP   = obsv_r(1:end-1,:);
    PP_p = obsv_r(1+1:end,:);
    threshold = 1e-7;
    PP_inv = pinv(PP,threshold);
    A_r = PP_inv*PP_p;
    
% D_r ~ !!!!! I'm assuming this is correct
    D_r = zeros(N_outputs,N_in);


%% Return the system
    sys = ss(A_r, B_r, C_r, D_r, SIM.Tsample);
end

%% %%%%%%%%%%%%%%%%%%%% OOOOOOOOOOOOOLLLDD %%%%%%%%%%%%%%%%%%%% 

%% Load Impulse Data
    % filename = getImpulseFilename(FLAG);
    % simsys = load(filename);

        % idx_vec = 1:simsys.SIM.TsMultiple:length(simsys.t_soln);
    % i_user_DT = simsys.i_user(idx_vec);

    % z_imp    = SIM.OutputMatrix * simsys.SV_soln';

    % t_imp = t_soln;

    % t_imp = simsys.t_soln;

    % save('F:\TylerFiles\~PhDWork\Presentations\Seminar\Spring2023\Pictures\HoKalmanROM\HKData.mat','simsys','z_imp','t_imp','g_k','idx_vec','pulse_idx','i_user_DT')


%% Plot Impulse Response
    % if FLAG.Analysis.PlotImp
    %     ODE_data_adj = inv(diag(range))*(z_imp -IC);
    %     for i = 1:N.DesOut
    %         figure
    %         hold on
    %         % plot(t_imp                            ,z_imp(i,:)/multiple(i,i) ,'-ok'                      ,'DisplayName','ODE')
    %         plot(t_imp                            ,ODE_data_adj(i,:) ,'-ok'                      ,'DisplayName','ODE')
    %         plot(t_imp(idx_vec(pulse_idx-1:end),:),g_k(i,:)  , 'or','MarkerFaceColor','r','DisplayName','Selected Data') %+IC(i,1)
    %         title(['Data Used in Ho-Kalman Reduction ' RESULTS.Labels.title{i}])
    %         lgn = legend;
    %         lgn.Location = 'best';
    %         xlim([0,10])
    %     end
    % end


%% Initialize Function Outputs
    % sys      = cell(1,N.DesOut+1);
    % singVals = nan(1, N.DesOut+1);


%% Optimal
% if FLAG.UseOptimal
%     [ROMError_filename] = getROMErrorFilename(FLAG);
%     if isfile(ROMError_filename)
%         optimal_r_obj = load(ROMError_filename,'optimal_r');
%     else
%         optimal_r_obj = 0;
%     end
% end


        % for OO = 1:N.DesOut


            % if FLAG.UseInput_r
            %     r = SIM.Input_r;
            % else
            %     r = 10;
            % end

            % singVals(OO) = S(r,r);
    
            % C_r
            % ind_multiple_vec = [multiple(1,1) , multiple(OO,OO)];
            % ind_multiple     = diag(ind_multiple_vec);

                % PP   = obsv_r(1:end-N_outputs,:);
                % PP_p = obsv_r(N_outputs+1:end,:);




            % if FLAG.PlotSingVal
            %     figure
            %     semilogy(1:1:length(diag(S)) , diag(S),'-k','Linewidth',2)
            %     xlim([0,length(diag(S))])
            %     ylabel('Singular Values')
            %     title([RESULTS.Labels.title{OO} ' Singular Values of \Sigma for Hankel'])
            % end

