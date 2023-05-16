%% 
function Kalman_plotFnc(RESULTS,N,SIM,FLAG)
%%
FLAG.PLOT.NoNoiseCompare            = 0;    
    FLAG.PLOT.ode        = 1;
    FLAG.PLOT.SS_CT      = 0;
    FLAG.PLOT.SS_DT      = 0;
    FLAG.PLOT.ROM_HoKal  = 1;

FLAG.PLOT.NoisyPlant                = 0;
    FLAG.PLOT.NoNoise           = 0;
    FLAG.PLOT.NoisyPlant_CV     = 0; % (Q & R) , CV:  cell voltage
    FLAG.PLOT.NoisyPlant_All    = 1; % (Q)     , All: All Desired Outputs
    FLAG.PLOT.NoisyPlantROM_CV  = 0; % (Q & R) , CV:  cell voltage
    FLAG.PLOT.NoisyPlantROM_All = 0; % (Q)     , All: All Desired Outputs

FLAG.PLOT.EstimatorZ                = 1;
    FLAG.PLOT.ASYZ              = 0;
    FLAG.PLOT.VARZ              = 0; 
        FLAG.PLOT.PlantComparZ = 1;
    FLAG.All3                   = 1; % Plant, Asy, and Var
    

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
    if isfield(RESULTS,'Slink_plant')
        % No Noise Plant
        if isfield(RESULTS.Slink_plant,'z') && FLAG.PLOT.NoNoise
            figure
            plot(RESULTS.Slink_plant.z.t_soln , RESULTS.Slink_plant.z.z_soln,'k','LineWidth',2,'DisplayName','Slink plant')
            title('Simulink No Noise Plant Cell Voltage')
            xlabel('Time [s]')
            ylabel(RESULTS.Labels.unit{1})
        end

        % Noisy Plant
        if isfield(RESULTS,'Slink_plant') && FLAG.PLOT.NoisyPlant_CV
            figure
            plot(RESULTS.Slink_plant.zN.t_soln , RESULTS.Slink_plant.zN.z_soln,'k','LineWidth',2,'DisplayName','Slink Noisy plant')
            title('Noisy Plant Plant Cell Voltage')
            xlabel('Time [s]')
            ylabel(RESULTS.Labels.unit{1})
        end

        % Noisy Plant All Outputs
        if isfield(RESULTS,'Slink_plant') && FLAG.PLOT.NoisyPlant_All
            for i = 1:N.DesOut
                figure
                hold on
                plot(RESULTS.Slink_plant.zNA.t_soln , RESULTS.Slink_plant.zNA.z_soln(i,:),'k','LineWidth',2,'DisplayName','Slink Noisy plant')
                lgn = legend;
                lgn.Location = 'best';
                title([RESULTS.Labels.title{i} ' Noisy Output'])
                xlabel('Time [s]')
                ylabel(RESULTS.Labels.unit{i})
            end
        end

        % Compare Noisy Plant and Noisy ROM
        if isfield(RESULTS.Slink_plant, 'zN') && isfield(RESULTS.Slink_plant, 'zNR') && FLAG.PLOT.NoisyPlantROM_CV
            figure
            hold on
            plot(RESULTS.Slink_plant.zN.t_soln , RESULTS.Slink_plant.zN.z_soln              ,'o','LineWidth',2,'DisplayName','Slink Noisy Plant')
            plot(RESULTS.Slink_plant.zNR.t_soln, RESULTS.Slink_plant.zNR.z_soln + z_init(1) ,'k','LineWidth',2,'DisplayName','Slink Noisy ROM')
            lgn = legend;
            lgn.Location = 'best';
            title(['Compare Noisy Plant and ROM ' RESULTS.Labels.title{1}])
            xlabel('Time [s]')
            ylabel(RESULTS.Labels.unit{1})
        end

        % Compare Noisy Plant and Noisy ROM All Outputs
        if isfield(RESULTS.Slink_plant, 'zNA') && isfield(RESULTS.Slink_plant, 'zNRA') && FLAG.PLOT.NoisyPlantROM_All
            for i = 1:N.DesOut % (Noise isn't applied to the measurements, just process noise)
                figure
                hold on
                plot(RESULTS.Slink_plant.zNA.t_soln , RESULTS.Slink_plant.zNA.z_soln(i,:)              ,'o','LineWidth',2,'DisplayName','Slink Noisy Plant')
                plot(RESULTS.Slink_plant.zNRA.t_soln ,RESULTS.Slink_plant.zNRA.z_soln(:,i) + z_init(i) ,'k','LineWidth',2,'DisplayName','Slink Noisy ROM')
                lgn = legend;
                lgn.Location = 'best';
                title(['Compare Noisy Plant and ROM ' RESULTS.Labels.title{i}])
                xlabel('Time [s]')
                ylabel(RESULTS.Labels.unit{i})
            end
        end
    end
end


%% Estimator Z
% Asymtptotic Estimator
if FLAG.PLOT.EstimatorZ && isfield(RESULTS,'EST')
    if FLAG.PLOT.ASYZ && isfield(RESULTS.EST,'ASY')
        if isfield(RESULTS.EST.ASY,'z_soln_ALL')
            if FLAG.EstimatorModel == 1 || FLAG.EST.SepHK == 0
                [r,~] = size(RESULTS.EST.ASY.z_soln_ALL);
                for OO = 1:r
                    figure
                    hold on
                    plot(RESULTS.EST.ASY.t_soln , RESULTS.EST.ASY.z_soln_ALL{1}(OO,:),'or','LineWidth',2,'DisplayName',['Est ASY'])
                    if FLAG.PLOT.PlantComparZ && isfield(RESULTS.EST,'PLANT')
                        if OO == 1
                            plot(RESULTS.EST.PLANT.t_soln , RESULTS.EST.PLANT.z_soln{1}(OO,:),'k','LineWidth',2,'DisplayName',['Est Plant'])
                        else
                            plot(RESULTS.EST.PLANT.t_soln , RESULTS.EST.PLANT.z_soln_ALL{1}(OO,:),'k','LineWidth',2,'DisplayName',['Est Plant'])
                        end
                        title([RESULTS.Labels.title{OO} ' Compare Plant and Asymptotic Estimator Outputs'])
                    else
                        title([RESULTS.Labels.title{OO} ' Asymptotic Estimator Outputs'])
                    end
                    xlabel('Time [s]')
                    ylabel(RESULTS.Labels.unit{OO})
                    lgn = legend;
                    lgn.Location = 'best';
                end
            else
                for OO = 1:length(RESULTS.EST.ASY.z_soln_ALL)
                    figure
                    hold on
                    plot(RESULTS.EST.ASY.t_soln , RESULTS.EST.ASY.z_soln_ALL{OO}(2,:),'or','LineWidth',2,'DisplayName',['Est ASY'])
                    if FLAG.PLOT.PlantComparZ && isfield(RESULTS.EST,'PLANT')
                        plot(RESULTS.EST.PLANT.t_soln , RESULTS.EST.PLANT.z_soln_ALL{OO}(2,:),'k','LineWidth',2,'DisplayName',['Est Plant'])
                        title([RESULTS.Labels.title{OO} ' Compare Plant and Asymptotic Estimator Outputs'])
                    else
                        title([RESULTS.Labels.title{OO} ' Asymptotic Estimator Outputs'])
                    end
                    xlabel('Time [s]')
                    ylabel(RESULTS.Labels.unit{OO})
                    lgn = legend;
                    lgn.Location = 'best';
                end
            end
        end
    end
end

% Variable Estimator
if FLAG.PLOT.EstimatorZ && isfield(RESULTS,'EST')
    if FLAG.PLOT.VARZ && isfield(RESULTS.EST,'VAR')
        if isfield(RESULTS.EST.VAR,'z_soln_ALL')
            if FLAG.EstimatorModel == 1 || FLAG.EST.SepHK == 0
                [r,~] = size(RESULTS.EST.VAR.z_soln_ALL);
                for OO = 1:r
                    figure
                    hold on
                    plot(RESULTS.EST.VAR.t_soln , RESULTS.EST.VAR.z_soln_ALL{1}(OO,:),'og','LineWidth',2,'DisplayName',['Est VAR'])
                    if FLAG.PLOT.PlantComparZ && isfield(RESULTS.EST,'PLANT')
                        if OO == 1
                            plot(RESULTS.EST.PLANT.t_soln , RESULTS.EST.PLANT.z_soln{1}(OO,:),'k','LineWidth',2,'DisplayName',['Est Plant'])
                        else
                            plot(RESULTS.EST.PLANT.t_soln , RESULTS.EST.PLANT.z_soln_ALL{1}(OO,:),'k','LineWidth',2,'DisplayName',['Est Plant'])
                        end
                        title([RESULTS.Labels.title{OO} ' Compare Plant and Variable Estimator Outputs'])
                    else
                        title([RESULTS.Labels.title{OO} ' Variable Estimator Outputs'])
                    end
                    xlabel('Time [s]')
                    ylabel(RESULTS.Labels.unit{OO})
                    lgn = legend;
                    lgn.Location = 'best';
                end
            else
                for OO = 1:length(RESULTS.EST.VAR.z_soln_ALL)
                    figure
                    hold on
                    plot(RESULTS.EST.VAR.t_soln , RESULTS.EST.VAR.z_soln_ALL{OO}(2,:),'og','LineWidth',2,'DisplayName',['Est VAR'])
                    if FLAG.PLOT.PlantComparZ && isfield(RESULTS.EST,'PLANT')
                        plot(RESULTS.EST.PLANT.t_soln , RESULTS.EST.PLANT.z_soln_ALL{OO}(2,:),'k','LineWidth',2,'DisplayName',['Est Plant'])
                        title([RESULTS.Labels.title{OO} ' Compare Plant and Variable Estimator Outputs'])
                    else
                        title([RESULTS.Labels.title{OO} ' Variable Estimator Outputs'])
                    end
                    xlabel('Time [s]')
                    ylabel(RESULTS.Labels.unit{OO})
                    lgn = legend;
                    lgn.Location = 'best';
                end
            end
        end
    end
end

% All 3
if FLAG.PLOT.EstimatorZ && isfield(RESULTS,'EST')
    if FLAG.All3 && isfield(RESULTS.EST.VAR,'z_soln_ALL') && isfield(RESULTS.EST,'PLANT') && isfield(RESULTS.EST.ASY,'z_soln_ALL')
        if FLAG.EstimatorModel == 1 || FLAG.EST.SepHK == 0
            [r,~] = size( RESULTS.EST.VAR.z_soln_ALL{1} );
            for OO = 1:r
                figure
                hold on
                plot(RESULTS.EST.ASY.t_soln , RESULTS.EST.ASY.z_soln_ALL{1}(OO,:),'or','LineWidth',2,'DisplayName',['Asymptotic'])
                plot(RESULTS.EST.VAR.t_soln , RESULTS.EST.VAR.z_soln_ALL{1}(OO,:),'og','LineWidth',2,'DisplayName',['Variable'])
%                 if OO == 1
%                     plot(RESULTS.EST.PLANT.t_soln , RESULTS.EST.PLANT.z_soln{1}(OO,:),'k','LineWidth',2,'DisplayName',['Plant'])
%                 else
                    plot(RESULTS.EST.PLANT.t_soln , RESULTS.EST.PLANT.z_soln_ALL{1}(OO,:),'k','LineWidth',2,'DisplayName',['Plant'])
%                 end
                title([RESULTS.Labels.title{OO} ' Compare Plant to Variable and Asymptotic Estimators'])
                xlabel('Time [s]')
                ylabel(RESULTS.Labels.unit{OO})
                lgn = legend;
                lgn.Location = 'best';
            end
        else
            for OO = 1:length(RESULTS.EST.VAR.z_soln_ALL)
                figure
                hold on
                plot(RESULTS.EST.ASY.t_soln ,     RESULTS.EST.ASY.z_soln_ALL{OO}(2,:),'or','LineWidth',2,'DisplayName',['Asymptotic'])
                plot(RESULTS.EST.VAR.t_soln ,     RESULTS.EST.VAR.z_soln_ALL{OO}(2,:),'og','LineWidth',2,'DisplayName',['Variable'])
                plot(RESULTS.EST.PLANT.t_soln , RESULTS.EST.PLANT.z_soln_ALL{OO}(2,:),'k' ,'LineWidth',2,'DisplayName',['Plant'])
                title([RESULTS.Labels.title{OO} ' Compare Plant to Variable and Asymptotic Estimators'])
                xlabel('Time [s]')
                ylabel(RESULTS.Labels.unit{OO})
                lgn = legend;
                lgn.Location = 'best';
            end
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