function DAE_Circuit_plotFnc(filename)
load(filename)

%% FLAGS
FLAG.PLOT.NoNoiseCompare = 0;
    FLAG.PLOT.ode3       = 1;
    FLAG.PLOT.ode5       = 0;
    FLAG.PLOT.SS_CT3     = 1;
    FLAG.PLOT.SS_DT3     = 0;
    FLAG.PLOT.ROM_Mlab3  = 0;
    FLAG.PLOT.ROM_HoKal  = 1;

FLAG.PLOT.NoisyPlant     = 1;

FLAG.PLOT.EstimatorX        = 0;
    FLAG.PLOT.ASYX          = 1;
    FLAG.PLOT.VARX          = 0;
    FLAG.PLOT.PlantComparX  = 0; % This is C_cross on plant outputs to guess what the states would be. and its not good
    FLAG.PLOT.ROMAllComparX = 1; % This is the ROM used in estimator but with a noisy input
    FLAG.PLOT.SlinkAllComparX = 0; % This is all 3 states

FLAG.PLOT.EstimatorZ       = 1;
    FLAG.PLOT.ASYZ         = 1;
    FLAG.PLOT.VARZ         = 0;
    FLAG.PLOT.PlantComparZ = 1;
    FLAG.PLOT.ROMAllComparZ = 0;
    FLAG.PLOT.SlinkAllComparZ = 1;

%FLAG.PLOT.ErrorCalcs        = 0;
FLAG.PLOT.GenComparData     = 0;
FLAG.PLOT.ComparSVD2Pinf    = 0;


%% No Noise Comparison
if FLAG.PLOT.NoNoiseCompare
    figure
    hold on
    if FLAG.PLOT.ROM_HoKal && isfield(RESULTS,'ROM_HoKal')
        plot(RESULTS.ROM_HoKal.t_soln,RESULTS.ROM_HoKal.z_soln(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1 ROM Ho-Kal')
        plot(RESULTS.ROM_HoKal.t_soln,RESULTS.ROM_HoKal.z_soln(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2 ROM Ho-Kal')
        plot(RESULTS.ROM_HoKal.t_soln,RESULTS.ROM_HoKal.z_soln(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3 ROM Ho-Kal')
    end
    if FLAG.PLOT.ROM_Mlab3 && isfield(RESULTS,'ROM_Mlab3')
        plot(RESULTS.ROM_Mlab3.t_soln,RESULTS.ROM_Mlab3.z_soln(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1 ROM Mlab3')
        plot(RESULTS.ROM_Mlab3.t_soln,RESULTS.ROM_Mlab3.z_soln(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2 ROM Mlab3')
        plot(RESULTS.ROM_Mlab3.t_soln,RESULTS.ROM_Mlab3.z_soln(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3 ROM Mlab3')
    end
    if FLAG.PLOT.SS_DT3    && isfield(RESULTS,'SS_DT3') 
        plot(RESULTS.SS_DT3.t_soln,RESULTS.SS_DT3.z_soln(:,P.V_1),'^','LineWidth',2,'DisplayName','V_1 SS DT3')
        plot(RESULTS.SS_DT3.t_soln,RESULTS.SS_DT3.z_soln(:,P.V_2),'^','LineWidth',2,'DisplayName','V_2 SS DT3')
        plot(RESULTS.SS_DT3.t_soln,RESULTS.SS_DT3.z_soln(:,P.V_3),'^','LineWidth',2,'DisplayName','V_3 SS DT3')
    end
    if FLAG.PLOT.SS_CT3    && isfield(RESULTS,'SS_CT3') 
        plot(RESULTS.SS_CT3.t_soln,RESULTS.SS_CT3.z_soln(:,P.V_1),'square','LineWidth',2,'DisplayName','V_1 SS CT3')
        plot(RESULTS.SS_CT3.t_soln,RESULTS.SS_CT3.z_soln(:,P.V_2),'square','LineWidth',2,'DisplayName','V_2 SS CT3')
        plot(RESULTS.SS_CT3.t_soln,RESULTS.SS_CT3.z_soln(:,P.V_3),'square','LineWidth',2,'DisplayName','V_3 SS CT3')
    end
    if FLAG.PLOT.ode3      &&  isfield(RESULTS,'ode3')
        plot(RESULTS.ode3.t_soln,RESULTS.ode3.x_soln(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1 ode3')
        plot(RESULTS.ode3.t_soln,RESULTS.ode3.x_soln(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2 ode3')
        plot(RESULTS.ode3.t_soln,RESULTS.ode3.x_soln(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3 ode3')
    end
    if FLAG.PLOT.ode5      &&  isfield(RESULTS,'ode5')
        plot(RESULTS.ode5.t_soln,RESULTS.ode5.x_soln(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1 ode5')
        plot(RESULTS.ode5.t_soln,RESULTS.ode5.x_soln(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2 ode5')
        plot(RESULTS.ode5.t_soln,RESULTS.ode5.x_soln(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3 ode5')
    end
    lgn = legend;
    title('No Noise Comparison')
    xlabel('Time [s]')
    ylabel('Voltage [V]')
end

%% Noisy Plant
if FLAG.PLOT.NoisyPlant
    % No Noise Plant
    if isfield(RESULTS,'Slink_plant')
        figure
        hold on
        [r,~] = size(RESULTS.Slink_plant.z.z_soln);
        for i = 1:r
            plot(RESULTS.Slink_plant.z.t_soln,RESULTS.Slink_plant.z.z_soln(i,:),'o','LineWidth',2,'DisplayName',['V_' num2str(i) ' Slink plant']) %%% Change name
        end
        title('Simulink No Noise Plant C_m')
        xlabel('Time [s]')
        ylabel('Voltage [V]')
    end
    
    % Noisy Plant
    if isfield(RESULTS,'Slink_plant')
        figure
        hold on
        [r,~] = size(RESULTS.Slink_plant.zN.z_soln);
        for i = 1:r
            plot(RESULTS.Slink_plant.zN.t_soln,RESULTS.Slink_plant.zN.z_soln(i,:),'k','LineWidth',2,'DisplayName',['V_' num2str(i) ' Slink Noisy plant'])
        end
        title('Noisy Plant C_m')
        xlabel('Time [s]')
        ylabel('Voltage [V]')
    end

    % Noisy Plant All Outputs
    if isfield(RESULTS,'Slink_plant')
        figure
        hold on
        [r,~] = size(RESULTS.Slink_plant.zNA.z_soln); %zNA is measurement Noisy All (Noise isn't applied to the measurements, just process noise)
        for i = 1:r
            plot(RESULTS.Slink_plant.zNA.t_soln,RESULTS.Slink_plant.zNA.z_soln(i,:),'k','LineWidth',2,'DisplayName',['V_' num2str(i) ' Slink Noisy plant'])
        end
        title('Noisy Plant C_{all}')
        xlabel('Time [s]')
        ylabel('Voltage [V]')
    end

    % Compare Noisy Plant and Noisy ROM
    if isfield(RESULTS,'Slink_plant')
        figure
        hold on
        % PLant
        [r,~] = size(RESULTS.Slink_plant.zNA.z_soln); %zNA is measurement Noisy All (Noise isn't applied to the measurements, just process noise)
        for i = 1:r
            plot(RESULTS.Slink_plant.zNA.t_soln,RESULTS.Slink_plant.zNA.z_soln(i,:),'o','LineWidth',2,'DisplayName',['V_' num2str(i) ' Slink Noisy plant'])
        end
        % ROM
        [r,~] = size(RESULTS.Slink_plant.zNRA.z_soln); %zNA is measurement Noisy All (Noise isn't applied to the measurements, just process noise)
        for i = 1:r
            plot(RESULTS.Slink_plant.zNRA.t_soln,RESULTS.Slink_plant.zNRA.z_soln(i,:),'k','LineWidth',2,'DisplayName',['V_' num2str(i) ' Slink Noisy plant'])
        end
        title('Compare Plant and ROM C_{all}')
        xlabel('Time [s]')
        ylabel('Voltage [V]')

    end
end

%% Estimator States
if FLAG.PLOT.EstimatorX && isfield(RESULTS,'EST')
    if FLAG.PLOT.ASYX && isfield(RESULTS.EST,'ASY')
        figure
        hold on
        [r,~] = size(RESULTS.EST.ASY.x_soln);
        for i = 1:r
            plot(RESULTS.EST.ASY.t_soln , RESULTS.EST.ASY.x_soln(i,:),'o','LineWidth',2,'DisplayName',['x_' num2str(i) ' Est ASY'])
        end
        if FLAG.PLOT.PlantComparX && isfield(RESULTS.EST,'PLANT')
            for i = 1:r
                plot(RESULTS.EST.PLANT.t_soln , RESULTS.EST.PLANT.x_soln(i,:),'k','LineWidth',2,'DisplayName',['x_' num2str(i) ' Est Plant'])
            end
            title('Compare Plant and Asymptotic Estimator States')
        else
            title('Asymptotic Estimator States')
        end
        if FLAG.PLOT.ROMAllComparX && isfield(RESULTS.Slink_plant,'xNR')
            for i = 1:r
                plot(RESULTS.Slink_plant.xNR.t_soln , RESULTS.Slink_plant.xNR.x_soln(i,:),'k','LineWidth',2,'DisplayName',['x_' num2str(i) ' Noisy ROM'])
            end
            %title('Compare Plant and Asymptotic Estimator States')
        else
            %title('Asymptotic Estimator States')
        end
        if FLAG.PLOT.SlinkAllComparX && isfield(RESULTS.Slink_plant,'xN')
            for i = 1:r
                plot(RESULTS.Slink_plant.xN.t_soln , RESULTS.Slink_plant.xN.x_soln(i,:),'or','LineWidth',2,'DisplayName',['x_' num2str(i) ' Noisy Plant'])
            end
            %title('Compare Plant and Asymptotic Estimator States')
        else
            %title('Asymptotic Estimator States')
        end
        xlabel('Time [s]')
        ylabel('Voltage [V]')
        lgn = legend;
    end


    if FLAG.PLOT.VARX && isfield(RESULTS.EST,'VAR')
        figure
        hold on
        [r,~] = size(RESULTS.EST.VAR.x_soln);
        for i = 1:r
            plot(RESULTS.EST.VAR.t_soln , RESULTS.EST.VAR.x_soln(i,:),'o','LineWidth',2,'DisplayName',['x_' num2str(i) ' Est VAR'])
        end
        if FLAG.PLOT.PlantComparX && isfield(RESULTS.EST,'PLANT')
            for i = 1:r
                plot(RESULTS.EST.PLANT.t_soln , RESULTS.EST.PLANT.x_soln(i,:),'k','LineWidth',2,'DisplayName',['x_' num2str(i) ' Est Plant'])
            end
            title('Compare Plant and Variable Estimator States')
        else
            title('Variable Estimator States')
        end
        if FLAG.PLOT.ROMAllComparX && isfield(RESULTS.Slink_plant,'xNR')
            for i = 1:r
                plot(RESULTS.Slink_plant.xNR.t_soln , RESULTS.Slink_plant.xNR.x_soln(i,:),'k','LineWidth',2,'DisplayName',['x_' num2str(i) ' Noisy ROM'])
            end
            %title('Compare Plant and Asymptotic Estimator States')
        else
            %title('Asymptotic Estimator States')
        end
        if FLAG.PLOT.SlinkAllComparX && isfield(RESULTS.Slink_plant,'xN')
            for i = 1:r
                plot(RESULTS.Slink_plant.xN.t_soln , RESULTS.Slink_plant.xN.x_soln(i,:),'or','LineWidth',2,'DisplayName',['x_' num2str(i) ' Noisy Plant'])
            end
            %title('Compare Plant and Asymptotic Estimator States')
        else
            %title('Asymptotic Estimator States')
        end
        xlabel('Time [s]')
        ylabel('Voltage [V]')
        lgn = legend;
    end
    
end

%% Estimator
if FLAG.PLOT.EstimatorZ && isfield(RESULTS,'EST')
    if FLAG.PLOT.ASYZ && isfield(RESULTS.EST,'ASY')
        figure
        hold on
        [r,~] = size(RESULTS.EST.VAR.z_soln);
        for i = 1:r
            plot(RESULTS.EST.ASY.t_soln , RESULTS.EST.ASY.z_soln(i,:),'o','LineWidth',2,'DisplayName',['V_' num2str(i) ' Est ASY'])
        end
        if FLAG.PLOT.PlantComparZ && isfield(RESULTS.EST,'PLANT')
            for i = 1:r
                plot(RESULTS.EST.PLANT.t_soln , RESULTS.EST.PLANT.z_soln(i,:),'k','LineWidth',2,'DisplayName',['V_' num2str(i) ' Est Plant'])
            end
            title('Compare Plant and Asymptotic Estimator Outputs')
        else
            title('Asymptotic Estimator Outputs')
        end
        if FLAG.PLOT.ROMAllComparZ && isfield(RESULTS.Slink_plant,'zNRA')
            for i = 1:r
                plot(RESULTS.Slink_plant.zNRA.t_soln , RESULTS.Slink_plant.zNRA.z_soln(i,:),'k','LineWidth',2,'DisplayName',['V_' num2str(i) ' Noisy ROM Output'])
            end
            %title('Compare Plant and Asymptotic Estimator States')
        else
            %title('Asymptotic Estimator States')
        end
        if FLAG.PLOT.SlinkAllComparZ && isfield(RESULTS.Slink_plant,'zNA')
            for i = 1:r
                plot(RESULTS.Slink_plant.zNA.t_soln , RESULTS.Slink_plant.zNA.z_soln(i,:),'r','LineWidth',2,'DisplayName',['V_' num2str(i) ' Noisy Plant Output'])
            end
            %title('Compare Plant and Asymptotic Estimator States')
        else
            %title('Asymptotic Estimator States')
        end
        xlabel('Time [s]')
        ylabel('Voltage [V]')
        lgn = legend;

        % Actually compare all outputs
        if FLAG.PLOT.SlinkAllComparZ && isfield(RESULTS.Slink_plant,'zNA')
            figure
            hold on
            [r,~] = size(RESULTS.Slink_plant.zNA.z_soln);
            for i = 1:r
                plot(RESULTS.Slink_plant.zNA.t_soln , RESULTS.Slink_plant.zNA.z_soln(i,:),'o','LineWidth',2,'DisplayName',['V_' num2str(i) ' Noisy Plant Output'])
            end
            [r,~] = size(RESULTS.EST.ASY.z_soln_ALL);
            for i = 1:r
                plot(RESULTS.EST.ASY.t_soln , RESULTS.EST.ASY.z_soln_ALL(i,:),'k','LineWidth',2,'DisplayName',['V_' num2str(i) ' Estimator Output'])
            end
            title('Compare All Outputs')
            xlabel('Time [s]')
            ylabel('Voltage [V]')
            lgn = legend;
        end
    end


    if FLAG.PLOT.VARZ && isfield(RESULTS.EST,'VAR')
        figure
        hold on
        [r,~] = size(RESULTS.EST.VAR.z_soln);
        for i = 1:r
            plot(RESULTS.EST.VAR.t_soln , RESULTS.EST.VAR.z_soln(i,:),'o','LineWidth',2,'DisplayName',['V_' num2str(i) ' Est VAR'])
        end
        if FLAG.PLOT.PlantComparZ && isfield(RESULTS.EST,'PLANT')
            for i = 1:r
                plot(RESULTS.EST.PLANT.t_soln , RESULTS.EST.PLANT.z_soln(i,:),'k','LineWidth',2,'DisplayName',['V_' num2str(i) ' Est Plant'])
            end
            title('Compare Plant and Variable Estimator Outputs')
        else
            title('Variable Estimator Outputs')
        end
        if FLAG.PLOT.ROMAllComparZ && isfield(RESULTS.Slink_plant,'zNRA')
            for i = 1:r
                plot(RESULTS.Slink_plant.zNRA.t_soln , RESULTS.Slink_plant.zNRA.z_soln(i,:),'k','LineWidth',2,'DisplayName',['V_' num2str(i) ' Noisy ROM Output'])
            end
            %title('Compare Plant and Asymptotic Estimator States')
        else
            %title('Asymptotic Estimator States')
        end
        if FLAG.PLOT.SlinkAllComparZ && isfield(RESULTS.Slink_plant,'zNA')
            for i = 1:r
                plot(RESULTS.Slink_plant.zNA.t_soln , RESULTS.Slink_plant.zNA.z_soln(i,:),'r','LineWidth',2,'DisplayName',['V_' num2str(i) ' Noisy Plant Output'])
            end
            %title('Compare Plant and Asymptotic Estimator States')
        else
            %title('Asymptotic Estimator States')
        end
        xlabel('Time [s]')
        ylabel('Voltage [V]')
        lgn = legend;
    end
end

%% Compare SVD and P_infty
if FLAG.PLOT.ComparSVD2Pinf

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

end