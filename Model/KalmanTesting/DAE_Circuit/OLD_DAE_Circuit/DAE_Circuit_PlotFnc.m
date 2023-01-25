%% DAE Circuit Plot Function
function DAE_Circuit_PlotFnc(filename)
%%
load(filename)


%% Plotting
% Voltage 1
    if FLAG.PLOT.V1Compare
        figure
        hold on
        title('States Step Response Voltage 1')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
    
        if FLAG.PLOT.DAE
            plot(t_DAE,x_DAE(:,P.V_1),'o','LineWidth',2,'DisplayName','DAE')
        end

        if FLAG.PLOT.Simulink_Noise
            plot(t_Slink_wN,x_Slink_wN(:,P.V_1),'LineWidth',2,'DisplayName','Plant wN')
        end

        if FLAG.PLOT.Simulink_NoNoise
            plot(t_Slink_NN,x_Slink_NN(:,P.V_1),'LineWidth',2,'DisplayName','Plant NN')
        end

        if FLAG.PLOT.Estimator
            if FLAG.PLOT.RunAsymtotic
                plot(t_asy,x_asy(:,P.V_1),'o','LineWidth',5,'DisplayName','Asymp')
            end
            if FLAG.PLOT.RunVariable
                plot(t_var,x_var(:,P.V_1),'o','LineWidth',2,'DisplayName','Variable K_k')
            end
            plot(t_asy,x_EST_sys(:,P.V_1),'LineWidth',2,'DisplayName','Plant')
        end
        lgn = legend;
        lgn.Location = 'southeast';
        xlim([0,t_final])
    end


%% Voltage 2
    if FLAG.PLOT.V2Compare
        figure
        hold on
        title('States Step Response Voltage 2')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
    
        if FLAG.PLOT.DAE
            plot(t_DAE,x_DAE(:,P.V_2),'o','LineWidth',2,'DisplayName','DAE')
        end

        if FLAG.PLOT.Simulink_Noise
            plot(t_Slink_wN,x_Slink_wN(:,P.V_2),'LineWidth',2,'DisplayName','Plant wN')
        end

        if FLAG.PLOT.Simulink_NoNoise
            plot(t_Slink_NN,x_Slink_NN(:,P.V_2),'LineWidth',2,'DisplayName','Plant NN')
        end

        if FLAG.PLOT.Estimator
            if FLAG.PLOT.RunAsymtotic
                plot(t_asy,x_asy(:,P.V_2),'o','LineWidth',5,'DisplayName','Asymp')
            end
            if FLAG.PLOT.RunVariable
                plot(t_var,x_var(:,P.V_2),'o','LineWidth',2,'DisplayName','Variable k')
            end
            plot(t_asy,x_EST_sys(:,P.V_2),'LineWidth',2,'DisplayName','Plant')
        end
        lgn = legend;
        xlim([0,t_final])
    end


%% Voltage 3
    if FLAG.PLOT.V3Compare
        figure
        hold on
        title('States Step Response Voltage 3')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        
        if FLAG.PLOT.DAE
            plot(t_DAE,x_DAE(:,P.V_3),'o','LineWidth',2,'DisplayName','DAE')
        end

        if FLAG.PLOT.Simulink_Noise
            plot(t_Slink_wN,x_Slink_wN(:,P.V_3),'LineWidth',2,'DisplayName','Plant wN')
        end

        if FLAG.PLOT.Simulink_NoNoise
            plot(t_Slink_NN,x_Slink_NN(:,P.V_3),'LineWidth',2,'DisplayName','Plant NN')
        end

        if FLAG.PLOT.Estimator
            if FLAG.PLOT.RunAsymtotic
                plot(t_asy,x_asy(:,P.V_3),'o','LineWidth',5,'DisplayName','Asymp')
            end
            if FLAG.PLOT.RunVariable
                plot(t_var,x_var(:,P.V_3),'o','LineWidth',2,'DisplayName','Variable k')
            end
            plot(t_asy,x_EST_sys(:,P.V_3),'LineWidth',2,'DisplayName','Plant')
        end
        lgn = legend;
        lgn.Location = 'southeast';
        xlim([0,t_final])
    end


%% Output Measurement Plotting
% Voltage 1
    if FLAG.PLOT.V1Compare && (FLAG.Measure == 4 || FLAG.Measure == 1)
        if FLAG.Measure == 1
            P.V_1 = 1;
        end        

        figure
        hold on
        title('Measurement Step Response Voltage 1')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
    
        if FLAG.PLOT.DAE
            plot(t_DAE,z_DAE(:,P.V_1),'o','LineWidth',2,'DisplayName','DAE')
        end
        if FLAG.PLOT.Simulink_Noise
            plot(t_Slink_wN,z_Slink_wN(:,P.V_1),'LineWidth',2,'DisplayName','Plant wN')
        end
        if FLAG.PLOT.Simulink_NoNoise
            plot(t_Slink_NN,z_Slink_NN(:,P.V_1),'LineWidth',2,'DisplayName','Plant NN')
        end

        if FLAG.PLOT.Estimator
            if FLAG.PLOT.RunAsymtotic
                plot(t_asy,z_asy(:,P.V_1),'o','LineWidth',5,'DisplayName','Asymp')
            end
            if FLAG.PLOT.RunVariable
                plot(t_var,z_var(:,P.V_1),'o','LineWidth',2,'DisplayName','Variable K_k')
            end
            plot(t_asy,z_EST_sys(:,P.V_1),'LineWidth',2,'DisplayName','Plant')
        end
        lgn = legend;
        lgn.Location = 'southeast';
        xlim([0,t_final])
    end


%% Voltage 2
    if FLAG.PLOT.V2Compare && (FLAG.Measure == 4 || FLAG.Measure == 2)
        if FLAG.Measure == 2
            P.V_2 = 1;
        end

        figure
        hold on
        title('Measurement Step Response Voltage 2')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
    
        if FLAG.PLOT.DAE
            plot(t_DAE,z_DAE(:,P.V_2),'o','LineWidth',2,'DisplayName','DAE')
        end
        if FLAG.PLOT.Simulink_Noise
            plot(t_Slink_wN,z_Slink_wN(:,P.V_2),'LineWidth',2,'DisplayName','Plant wN')
        end
        if FLAG.PLOT.Simulink_NoNoise
            plot(t_Slink_NN,z_Slink_NN(:,P.V_2),'LineWidth',2,'DisplayName','Plant NN')
        end

        if FLAG.PLOT.Estimator
            if FLAG.PLOT.RunAsymtotic
                plot(t_asy,z_asy(:,P.V_2),'o','LineWidth',5,'DisplayName','Asymp')
            end
            if FLAG.PLOT.RunVariable
                plot(t_var,z_var(:,P.V_2),'o','LineWidth',2,'DisplayName','Variable k')
            end
            plot(t_asy,z_EST_sys(:,P.V_2),'LineWidth',2,'DisplayName','Plant')
        end
        lgn = legend;
        xlim([0,t_final])
    end

    
%% Voltage 3
    if FLAG.PLOT.V3Compare && (FLAG.Measure == 4 || FLAG.Measure == 3)
        if FLAG.Measure == 3
            P.V_3 = 1;
        end
        figure
        hold on
        title('Measurement Step Response Voltage 3')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        
        if FLAG.PLOT.DAE
            plot(t_DAE,z_DAE(:,P.V_3),'o','LineWidth',2,'DisplayName','DAE')
        end
        if FLAG.PLOT.Simulink_Noise
            plot(t_Slink_wN,z_Slink_wN(:,P.V_3),'LineWidth',2,'DisplayName','Plant wN')
        end
        if FLAG.PLOT.Simulink_NoNoise
            plot(t_Slink_NN,z_Slink_NN(:,P.V_3),'LineWidth',2,'DisplayName','Plant NN')
        end

        if FLAG.PLOT.Estimator
            if FLAG.PLOT.RunAsymtotic
                plot(t_asy,z_asy(:,P.V_3),'o','LineWidth',5,'DisplayName','Asymp')
            end
            if FLAG.PLOT.RunVariable
                plot(t_var,z_var(:,P.V_3),'o','LineWidth',2,'DisplayName','Variable k')
            end
            plot(t_asy,z_EST_sys(:,P.V_3),'LineWidth',2,'DisplayName','Plant')
        end
        lgn = legend;
        lgn.Location = 'southeast';
        xlim([0,t_final])
    end


%% Arrange Figures
FigArrange = 1;
fig = gcf;
NumFig = fig.Number;
if FigArrange == 1
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
