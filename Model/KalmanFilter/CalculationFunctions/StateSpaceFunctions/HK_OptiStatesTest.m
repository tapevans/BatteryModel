function HK_OptiStatesTest(SIM,N,P,FLAG,RESULTS)
%% Load ODE Data
% Tswitch10
FLAG.InputMode = 5;
FLAG.PRBSAmp   = 1;
FLAG.Tswitch   = 10;
[Tswitch10.t_soln,  Tswitch10.x_soln,  Tswitch10.z_soln,  i_user_time, i_user_value] = getNoNoiseODE(SIM,FLAG,N,P);
% Input Signal
t_vec  = (0 : SIM.Ts : i_user_time(end))';
i_user = interp1(i_user_time , i_user_value , t_vec);
% Tswitch10.NewInputSignal = [t_vec , i_user];
Tswitch10.NewInputSignal = SIM.InputSignal;

if 0
    figure
    hold on
    plot(i_user_time        ,i_user_value                               ,'k','LineWidth' ,2,'DisplayName','ODE')
    plot(Tswitch10.NewInputSignal(:,1)  ,Tswitch10.NewInputSignal(:,2)  ,'ro','LineWidth',2,'DisplayName','DT')
    xlabel('Time [s]')
    ylabel('Current [A m^{-2}]')
    title('Compare ODE to DT Input Signal')
    lgn = legend;
end

OO = 7;
RR_end = 20;

for RR = 1:RR_end
    %% Create ROM
    SIM.Input_r = RR;
    [sys_HK,singVals] = getHoKalmanROM(SIM,N,P,FLAG,RESULTS);


    %% Get Singular Values
    %     for OO = 1:N.DesOut+1
    %         singularValue(1,RR,OO) = singVals(1,OO);
    %     end


    %% Get ROM simulation data for each test for all desired outputs
    %     for OO = 1:N.DesOut+1
    %% Step Data
    %         if OO < N.DesOut+1
    %             x_red = zeros(length(sys_HK{OO}.A),1);
    %             [z_soln, t_soln, ~] = lsim(sys_HK{OO},step.NewInputSignal(:,2),step.NewInputSignal(:,1),x_red);
    %             if OO == 1
    %                 step.individualOutputs.t_soln   = t_soln;
    %             end
    %             step.individualOutputs.z_soln(:,OO) = z_soln(:,2);
    %         else
    %             x_red = zeros(length(sys_HK{end}.A),1);
    %             [step.combinedOutputs.z_soln, step.combinedOutputs.t_soln, step.combinedOutputs.x_soln] = lsim(sys_HK{end},step.NewInputSignal(:,2),step.NewInputSignal(:,1),x_red);
    %         end

    %% Tswitch10
    if OO < N.DesOut+1
        x_red = zeros(length(sys_HK{OO}.A),1);
        [z_soln, t_soln, ~] = lsim(sys_HK{OO} , Tswitch10.NewInputSignal(:,2) , Tswitch10.NewInputSignal(:,1) , x_red);
        % if OO == 1
        Tswitch10.individualOutputs.t_soln = t_soln;
        % end
        Tswitch10.individualOutputs.z_soln = z_soln;
    else
        x_red = zeros(length(sys_HK{end}.A),1);
        [Tswitch10.combinedOutputs.z_soln, Tswitch10.combinedOutputs.t_soln, Tswitch10.combinedOutputs.x_soln] = lsim(sys_HK{end},Tswitch10.NewInputSignal(:,2),Tswitch10.NewInputSignal(:,1),x_red);
    end
    % x_red = zeros(length(sys_HK{OO}.A),1);
    % [z_soln, t_soln, ~] = lsim(sys_HK{OO},Tswitch10.NewInputSignal(:,2),Tswitch10.NewInputSignal(:,1),x_red);
    % Tswitch10.individualOutputs.t_soln = t_soln;
    % Tswitch10.individualOutputs.z_soln = z_soln;

    %% Tswitch100
    %         if OO < N.DesOut+1
    %             x_red = zeros(length(sys_HK{OO}.A),1);
    %             [z_soln, t_soln, ~] = lsim(sys_HK{OO},Tswitch100.NewInputSignal(:,2),Tswitch100.NewInputSignal(:,1),x_red);
    %             if OO == 1
    %                 Tswitch100.individualOutputs.t_soln   = t_soln;
    %             end
    %             Tswitch100.individualOutputs.z_soln(:,OO) = z_soln(:,2);
    %         else
    %             x_red = zeros(length(sys_HK{end}.A),1);
    %             [Tswitch100.combinedOutputs.z_soln, Tswitch100.combinedOutputs.t_soln, Tswitch100.combinedOutputs.x_soln] = lsim(sys_HK{end},Tswitch100.NewInputSignal(:,2),Tswitch100.NewInputSignal(:,1),x_red);
    %         end
    %     end
    %% Add Initial Offset back to SS
    %     step.individualOutputs.z_soln = step.individualOutputs.z_soln + SIM.y_0_FOM';
    %     step.combinedOutputs.z_soln   = step.combinedOutputs.z_soln   + SIM.y_0_FOM';

    if OO ~= 7
        Tswitch10.individualOutputs.z_soln = Tswitch10.individualOutputs.z_soln + SIM.y_0_FOM([1,OO],1)';
    else
        Tswitch10.combinedOutputs.z_soln   = Tswitch10.combinedOutputs.z_soln + SIM.y_0_FOM';
    end
    %     Tswitch10.combinedOutputs.z_soln   = Tswitch10.combinedOutputs.z_soln   + SIM.y_0_FOM';

    %     Tswitch100.individualOutputs.z_soln = Tswitch100.individualOutputs.z_soln + SIM.y_0_FOM';
    %     Tswitch100.combinedOutputs.z_soln   = Tswitch100.combinedOutputs.z_soln   + SIM.y_0_FOM';


    %% Shift Data
    % Tswitch10.individualOutputs.t_soln = Tswitch10.individualOutputs.t_soln(1:end-2);
    % Tswitch10.individualOutputs.z_soln = Tswitch10.individualOutputs.z_soln(3:end,:);

    %% Calculate Interpolated ODE Data
    for OOO = 1:N.DesOut
        %     % Step Data
        %         step.z_interpODE(OO,:) = interp1(step.t_soln, step.z_soln(OO,:) ,step.NewInputSignal(:,1));

        % Tswitch10
        Tswitch10.z_interpODE(OOO,:) = interp1(Tswitch10.t_soln, Tswitch10.z_soln(OOO,:) ,Tswitch10.NewInputSignal(:,1));

        %     % Tswitch100
        %         Tswitch100.z_interpODE(OO,:) = interp1(Tswitch100.t_soln, Tswitch100.z_soln(OO,:) ,Tswitch100.NewInputSignal(:,1));
    end

    if 0
        if RR >=5 && RR < 11
            if OO < N.DesOut+1
                figure
                hold on
                plot(Tswitch10.t_soln               ,Tswitch10.z_soln(1,:)                    , 'k','LineWidth',2,'DisplayName','ODE')
                plot(Tswitch10.NewInputSignal(:,1)  ,Tswitch10.z_interpODE(1,:)               ,'bo','LineWidth',2,'DisplayName','ODE interp')
                plot(Tswitch10.NewInputSignal(:,1)  ,Tswitch10.individualOutputs.z_soln(:,1)  ,'ro','LineWidth',5,'DisplayName','DT')
                xlabel('Time [s]')
                ylabel('Voltage [V]')
                title(['Compare ODE to DT , r= ' num2str(RR)])
                lgn = legend;
                xlim([50,80])

                figure
                hold on
                plot(Tswitch10.t_soln               ,Tswitch10.z_soln(OO,:)                   , 'k','LineWidth',2,'DisplayName','ODE')
                plot(Tswitch10.NewInputSignal(:,1)  ,Tswitch10.z_interpODE(OO,:)              ,'bo','LineWidth',2,'DisplayName','ODE interp')
                plot(Tswitch10.NewInputSignal(:,1)  ,Tswitch10.individualOutputs.z_soln(:,2)  ,'ro','LineWidth',5,'DisplayName','DT')
                xlabel('Time [s]')
                ylabel('Other')
                title(['Compare ODE to DT , r= ' num2str(RR)])
                lgn = legend;
                xlim([50,80])
            else
                for OOO = 1:N.DesOut
                    figure
                    hold on
                    plot(Tswitch10.t_soln               ,Tswitch10.z_soln(OOO,:)                    , 'k','LineWidth',2,'DisplayName','ODE')
                    plot(Tswitch10.NewInputSignal(:,1)  ,Tswitch10.z_interpODE(OOO,:)               ,'bo','LineWidth',2,'DisplayName','ODE interp')
                    plot(Tswitch10.NewInputSignal(:,1)  ,Tswitch10.combinedOutputs.z_soln(:,OOO)    ,'ro','LineWidth',2,'MarkerSize',10,'DisplayName','DT')
                    xlabel('Time [s]')
                    ylabel(RESULTS.Labels.unit{OOO})
                    title([RESULTS.Labels.title{OOO} ' Compare ODE to DT , r= ' num2str(RR) ])
                    lgn = legend;
                    xlim([430,470])
                end
            end
        end
    end


    %% Calculate Error at each time step
    % % Step Individual
    %     step.individualOutputs.Error       = abs(step.individualOutputs.z_soln'       - step.z_interpODE(:,1:end-2));
    % % Step Combined
    %     step.combinedOutputs.Error         = abs(step.combinedOutputs.z_soln'         - step.z_interpODE(:,1:end-2));
    if OO < N.DesOut+1
        % Tswitch10 Individual
        Tswitch10.individualOutputs.Error  = Tswitch10.individualOutputs.z_soln'  - Tswitch10.z_interpODE([1,OO],:);
    else
        % Tswitch10 Combined
        Tswitch10.combinedOutputs.Error    = Tswitch10.combinedOutputs.z_soln'    - Tswitch10.z_interpODE;
    end
    % % Tswitch100 Individual
    %     Tswitch100.individualOutputs.Error = abs(Tswitch100.individualOutputs.z_soln' - Tswitch100.z_interpODE(:,1:end-2));
    % % Tswitch100 Combined
    %     Tswitch100.combinedOutputs.Error   = abs(Tswitch100.combinedOutputs.z_soln'   - Tswitch100.z_interpODE(:,1:end-2));

    if OO < N.DesOut+1
        Tswitch10.individualOutputs.ErrorSQ    =    (Tswitch10.individualOutputs.Error).^2;
        Tswitch10.individualOutputs.ErrorSQSum = sum(Tswitch10.individualOutputs.ErrorSQ,2);
        Tswitch10.ind.saveErrorSqSum(:,RR)         = Tswitch10.individualOutputs.ErrorSQSum;
    else
        Tswitch10.combinedOutputs.ErrorSQ    =    (Tswitch10.combinedOutputs.Error).^2;
        Tswitch10.combinedOutputs.ErrorSQSum = sum(Tswitch10.combinedOutputs.ErrorSQ,2);
        Tswitch10.com.saveErrorSqSum(:,RR)   =     Tswitch10.combinedOutputs.ErrorSQSum;
    end


    %% Plots of the responses at different ranks
    if OO < N.DesOut+1
        if RR >=5 && RR < 11
            figure
            title(['Rank: ' num2str(RR) 'VoltageLSQ:' num2str(Tswitch10.individualOutputs.ErrorSQSum(1)) 'OtherLSQ:' num2str(Tswitch10.individualOutputs.ErrorSQSum(2))])

            yyaxis left
            hold on
            plot(Tswitch10.individualOutputs.t_soln , Tswitch10.individualOutputs.z_soln(:,1),'o','Linewidth',2)
            % plot(Tswitch10.t_soln,Tswitch10.z_soln(1,:),'-k')
            plot(Tswitch10.NewInputSignal(:,1),Tswitch10.z_interpODE(1,:),'-k')
            ylabel('Voltage')
            hold off

            yyaxis right
            hold on
            plot(Tswitch10.individualOutputs.t_soln , Tswitch10.individualOutputs.z_soln(:,2),'o','Linewidth',2)
            % plot(Tswitch10.t_soln,Tswitch10.z_soln(OO,:),'-k')
            plot(Tswitch10.NewInputSignal(:,1),Tswitch10.z_interpODE(OO,:),'-k')
            ylabel('Other')
            hold off
        end
    end

end


%% Plots of the error vs rank
if OO < N.DesOut+1
    [min_v , idx_volt] = min(Tswitch10.ind.saveErrorSqSum(1,:));
    [min_o , idx_othe] = min(Tswitch10.ind.saveErrorSqSum(2,:));
else
    [min_com , idx_com] = min(Tswitch10.com.saveErrorSqSum,[],2);
end

if OO < N.DesOut+1
    figure
    % plot(1:1:length(Tswitch10.ind.saveErrorSqSum) , Tswitch10.ind.saveErrorSqSum(1,:),'-k','Linewidth',2)
    semilogy(1:1:length(Tswitch10.ind.saveErrorSqSum) , Tswitch10.ind.saveErrorSqSum(1,:),'-k','Linewidth',2)
    hold on
    plot(idx_volt,min_v,'ro','Linewidth',2,'MarkerSize',10,'DisplayName','Voltage Minimum')
    plot(idx_othe,Tswitch10.ind.saveErrorSqSum(1,idx_othe),'bo','Linewidth',2,'MarkerSize',10,'DisplayName',[RESULTS.Labels.title{OO} ' Minimum'])
    xlabel('Number of States')
    ylabel(RESULTS.Labels.unit{1})
    title([RESULTS.Labels.title{1} ' \Sigma error^2'])
    xlim([1,RR_end])
    % ylim([0 1e-3])
    lgn = legend;

    figure
    % plot(1:1:length(Tswitch10.ind.saveErrorSqSum) , Tswitch10.ind.saveErrorSqSum(2,:),'-k','Linewidth',2)
    semilogy(1:1:length(Tswitch10.ind.saveErrorSqSum) , Tswitch10.ind.saveErrorSqSum(2,:),'-k','Linewidth',2)
    hold on
    plot(idx_othe,min_o,'bo','Linewidth',2,'MarkerSize',10)
    plot(idx_volt,Tswitch10.ind.saveErrorSqSum(2,idx_volt),'ro','Linewidth',2,'MarkerSize',10)
    xlabel('Number of States')
    ylabel(RESULTS.Labels.unit{OO})
    title([RESULTS.Labels.title{OO} ' \Sigma error^2'])
    xlim([1,RR_end])
    % ylim([0 1e-1])
else
    for i = 1:N.DesOut
        figure
        % plot(1:1:RR_end , Tswitch10.com.saveErrorSqSum(1,:),'-k','Linewidth',2)
        semilogy(1:1:RR_end , Tswitch10.com.saveErrorSqSum(i,:),'-k','Linewidth',2)
        hold on
        % plot(idx_volt,min_v,'ro','Linewidth',2,'MarkerSize',10,'DisplayName','Voltage Minimum')
        plot(idx_com(i),min_com(i),'bo','Linewidth',2,'MarkerSize',10,'DisplayName',[RESULTS.Labels.title{i} ' Minimum @ ' num2str(idx_com(i))])
        xlabel('Number of States')
        ylabel(RESULTS.Labels.unit{i})
        title([RESULTS.Labels.title{i} ' \Sigma error^2' ', Minimum @ ' num2str(idx_com(i))])
        xlim([1,RR_end])
        % ylim([0 1e-3])
        % lgn = legend;
    end
end


%% Arrange Figures
FigArrange = 1;
if FigArrange == 1
    fig = gcf;
    NumFig = fig.Number;

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


end % Function