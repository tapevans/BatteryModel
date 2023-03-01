%% 
function Kalman_plotFnc(RESULTS,N,SIM)
%%
FLAG.PLOT.NoNoiseCompare            = 1;    
    FLAG.PLOT.ode        = 1;
    FLAG.PLOT.SS_CT      = 0;
    FLAG.PLOT.SS_DT      = 0;
    FLAG.PLOT.ROM_HoKal  = 1;

FLAG.PLOT.NoisyPlant                = 0;
    FLAG.PLOT.NoNoise           = 0;
    FLAG.PLOT.NoisyPlant_CV     = 0; % (Q & R)
    FLAG.PLOT.NoisyPlant_All    = 0; % (Q)
    FLAG.PLOT.NoisyPlantROM_CV  = 0; % (Q & R)
    FLAG.PLOT.NoisyPlantROM_All = 1; % (Q)

FLAG.PLOT.EstimatorZ       = 1;
    FLAG.PLOT.ASYZ         = 1;
    FLAG.PLOT.VARZ         = 0; % Not implemented yet
    FLAG.PLOT.PlantComparZ = 1;

z_init = SIM.y_0_FOM;


%% No Noise Plots
if FLAG.PLOT.NoNoiseCompare
    if isfield(RESULTS,'ROM_HoKal') || isfield(RESULTS,'SS_DT') || isfield(RESULTS,'SS_CT') || isfield(RESULTS,'ode')
        for i = 1:N.DesOut
            figure
            hold on
            if FLAG.PLOT.ROM_HoKal && isfield(RESULTS,'ROM_HoKal')
                plot(RESULTS.ROM_HoKal.t_soln , RESULTS.ROM_HoKal.z_soln(:,i),'ro'     ,'LineWidth',2,'DisplayName','ROM Ho-Kal')
            end
            if FLAG.PLOT.SS_DT    && isfield(RESULTS,'SS_DT')
                plot(RESULTS.SS_DT.t_soln     , RESULTS.SS_DT.z_soln(:,i)    ,'b^'     ,'LineWidth',2,'DisplayName','SS DT')
            end
            if FLAG.PLOT.SS_CT    && isfield(RESULTS,'SS_CT')
                plot(RESULTS.SS_CT.t_soln     , RESULTS.SS_CT.z_soln(:,i)    ,'gsquare','LineWidth',2,'DisplayName','SS CT')
            end
            if FLAG.PLOT.ode      &&  isfield(RESULTS,'ode')
                plot(RESULTS.ode.t_soln       , RESULTS.ode.z_soln(i,:)      ,'k'      ,'LineWidth',2,'DisplayName','ODE')
            end
            lgn = legend;
            lgn.Location = 'best';
            title([RESULTS.Labels.title{i} ' No Noise Comparison'])
            xlabel('Time [s]')
            ylabel(RESULTS.Labels.unit{i})
    
        end
    end
end


%% Noisy Plant
if FLAG.PLOT.NoisyPlant
    % No Noise Plant
    if isfield(RESULTS,'Slink_plant') && FLAG.PLOT.NoNoise
        figure
        plot(RESULTS.Slink_plant.z.t_soln,RESULTS.Slink_plant.z.z_soln,'k','LineWidth',2,'DisplayName','Slink plant') 
        title('Simulink No Noise Plant Cell Voltage')
        xlabel('Time [s]')
        ylabel(RESULTS.Labels.unit{1})
    end
    
    % Noisy Plant
    if isfield(RESULTS,'Slink_plant') && FLAG.PLOT.NoisyPlant_CV
        figure
        plot(RESULTS.Slink_plant.zN.t_soln,RESULTS.Slink_plant.zN.z_soln,'k','LineWidth',2,'DisplayName','Slink Noisy plant')
        title('Noisy Plant Plant Cell Voltage')
        xlabel('Time [s]')
        ylabel(RESULTS.Labels.unit{1})
    end

    % Noisy Plant All Outputs
    if isfield(RESULTS,'Slink_plant') && FLAG.PLOT.NoisyPlant_All
        for i = 1:N.DesOut
            figure
            hold on
            plot(RESULTS.Slink_plant.zNA.t_soln,RESULTS.Slink_plant.zNA.z_soln(i,:),'k','LineWidth',2,'DisplayName','Slink Noisy plant')
            lgn = legend;
            lgn.Location = 'best';
            title([RESULTS.Labels.title{i} ' Noisy Output'])
            xlabel('Time [s]')
            ylabel(RESULTS.Labels.unit{i})
        end
    end

    % Compare Noisy Plant and Noisy ROM 
    if isfield(RESULTS,'Slink_plant') && FLAG.PLOT.NoisyPlantROM_CV
        figure
        hold on
        plot(RESULTS.Slink_plant.zN.t_soln ,RESULTS.Slink_plant.zN.z_soln              ,'o','LineWidth',2,'DisplayName','Slink Noisy Plant')
        plot(RESULTS.Slink_plant.zNR.t_soln,RESULTS.Slink_plant.zNR.z_soln + z_init(1) ,'k','LineWidth',2,'DisplayName','Slink Noisy ROM')
        lgn = legend;
        lgn.Location = 'best';
        title(['Compare Noisy Plant and ROM ' RESULTS.Labels.title{1}])
        xlabel('Time [s]')
        ylabel(RESULTS.Labels.unit{1})
    end

    % Compare Noisy Plant and Noisy ROM All Outputs
    if isfield(RESULTS,'Slink_plant') && FLAG.PLOT.NoisyPlantROM_All
        for i = 1:N.DesOut % (Noise isn't applied to the measurements, just process noise)
            figure
            hold on
            plot(RESULTS.Slink_plant.zNA.t_soln ,RESULTS.Slink_plant.zNA.z_soln(i,:)              ,'o','LineWidth',2,'DisplayName','Slink Noisy Plant')
            plot(RESULTS.Slink_plant.zNRA.t_soln,RESULTS.Slink_plant.zNRA.z_soln(:,i) + z_init(i) ,'k','LineWidth',2,'DisplayName','Slink Noisy ROM')
            lgn = legend;
            lgn.Location = 'best';
            title(['Compare Noisy Plant and ROM ' RESULTS.Labels.title{i}])
            xlabel('Time [s]')
            ylabel(RESULTS.Labels.unit{i})
        end
    end
end


%% Estimator Z
if FLAG.PLOT.EstimatorZ && isfield(RESULTS,'EST')
    if FLAG.PLOT.ASYZ && isfield(RESULTS.EST,'ASY')
        [r,~] = size(RESULTS.EST.ASY.z_soln_ALL);
        for i = 1:r
            figure
            hold on
            plot(RESULTS.EST.ASY.t_soln , RESULTS.EST.ASY.z_soln_ALL(i,:),'o','LineWidth',2,'DisplayName',['Est ASY'])
            if FLAG.PLOT.PlantComparZ && isfield(RESULTS.EST,'PLANT')
                plot(RESULTS.EST.PLANT.t_soln , RESULTS.EST.PLANT.z_soln_ALL(i,:),'k','LineWidth',2,'DisplayName',['Est Plant'])
                title([RESULTS.Labels.title{i} ' Compare Plant and Asymptotic Estimator Outputs'])
            else
                title([RESULTS.Labels.title{i} ' Asymptotic Estimator Outputs'])
            end
            xlabel('Time [s]')
            ylabel(RESULTS.Labels.unit{i})
        end
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
end