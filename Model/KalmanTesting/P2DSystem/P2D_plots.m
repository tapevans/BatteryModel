%% Plot the results of P2D

FLAG.PlotResults = 1;

FLAG.Plot.CellVoltage = 1;
FLAG.Plot.DelPhi      = 1;
FLAG.Plot.C_Liion     = 1;
FLAG.Plot.X_surf      = 1;
FLAG.Plot.i_Far       = 1;
FLAG.Plot.eta         = 1;

FLAG.Plot.sim_impl = 0;
FLAG.Plot.sim_sine = 0;
FLAG.Plot.ROM_step = 0;
FLAG.Plot.sim_step = 0;
FLAG.Plot.Slink.plant       = 0;
FLAG.Plot.Slink.plant_noise = 1; 
FLAG.Plot.Slink.ROM         = 0;
FLAG.Plot.Slink.Est         = 1;
FLAG.Plot.Slink.Est_asy     = 0;

%%
cell_of_interest = sim_data.step.N.CV_Region_AN(end);
if FLAG.Plot.sim_impl
    xmin  = 0;
    limit = 5;
else
    xmin  = 0;
    limit = Inf;

%     xmin = 2-1e-3;
%     limit = 1.2e-4;
%     limit = 20;
end

start_idx_est     = 1;
start_idx_est_asy = 2;


%%
if FLAG.PlotResults
    if FLAG.Plot.CellVoltage
        figure
        hold on
        title('Cell Voltage')
        if FLAG.Plot.Slink.plant_noise
            plot(t_slink_plant_noise, y_slink_plant_noise(:,localP.cell_voltage),'Linewidth',2,'DisplayName','Slink Plant Noise')
        end
        if FLAG.Plot.Slink.ROM
            plot(t_slink_ROM        , y_slink_ROM(:,localP.cell_voltage),'o','Linewidth',2,'DisplayName','Slink ROM')
        end
        if FLAG.Plot.Slink.Est
            plot(t_slink_est(start_idx_est:end), y_slink_est(start_idx_est:end,localP.cell_voltage),'o','Linewidth',2,'DisplayName','Slink Est')
        end
        if FLAG.Plot.Slink.Est_asy
            plot(t_slink_est_asy(start_idx_est_asy:end), y_slink_est_asy(start_idx_est_asy:end,localP.cell_voltage),'o','Linewidth',2,'DisplayName','Slink Est Asy')
        end        
        if FLAG.Plot.Slink.plant
            plot(t_slink_plant      , y_slink_plant(:,localP.cell_voltage),'Linewidth',2,'DisplayName','Slink Plant')
        end
        if FLAG.Plot.sim_impl
            plot(sim_data.impl.t_soln, sim_data.impl.cell_voltage,'Linewidth',2,'DisplayName','Impluse')
        end
        if FLAG.Plot.sim_sine
            plot(sim_data.sine.t_soln, sim_data.sine.cell_voltage,'Linewidth',2,'DisplayName','Harmonic')
        end
        if FLAG.Plot.ROM_step
            plot(t_SS_ROM,y_SS_ROM(:,localP.cell_voltage)+IC(localP.cell_voltage),'o','Linewidth',2,'DisplayName','ROM Step')
        end
        if FLAG.Plot.sim_step
            plot(sim_data.step.t_soln, sim_data.step.cell_voltage,'Linewidth',2,'DisplayName','Step')
        end
        
        lgn = legend;
        xlabel('Time [s]')
        ylabel('Voltage [V]')
        xlim([xmin,limit])
    end
    if FLAG.Plot.DelPhi
        figure
        hold on
        title('\Delta \phi')
        if FLAG.Plot.Slink.plant_noise
            plot(t_slink_plant_noise, y_slink_plant_noise(:,localP.delta_phi),'Linewidth',2,'DisplayName','Slink Plant Noise')
        end
        if FLAG.Plot.Slink.ROM
            plot(t_slink_ROM        , y_slink_ROM(:,localP.delta_phi),'o','Linewidth',2,'DisplayName','Slink ROM')
        end
        if FLAG.Plot.Slink.Est
            plot(t_slink_est(start_idx_est:end), y_slink_est(start_idx_est:end,localP.delta_phi),'o','Linewidth',2,'DisplayName','Slink Est')
        end
        if FLAG.Plot.Slink.Est_asy
            plot(t_slink_est_asy(start_idx_est_asy:end), y_slink_est_asy(start_idx_est_asy:end,localP.delta_phi),'o','Linewidth',2,'DisplayName','Slink Est Asy')
        end        
        if FLAG.Plot.Slink.plant
            plot(t_slink_plant      , y_slink_plant(:,localP.delta_phi),'Linewidth',2,'DisplayName','Slink Plant')
        end
        if FLAG.Plot.sim_impl
            plot(sim_data.impl.t_soln, sim_data.impl.del_phi(:,cell_of_interest),'Linewidth',2,'DisplayName','Impluse')
        end
        if FLAG.Plot.sim_sine
            plot(sim_data.sine.t_soln, sim_data.sine.del_phi(:,cell_of_interest),'Linewidth',2,'DisplayName','Harmonic')
        end
        if FLAG.Plot.ROM_step
            plot(t_SS_ROM,y_SS_ROM(:,localP.delta_phi)+IC(localP.delta_phi),'o','Linewidth',2,'DisplayName','ROM Step')
        end
        if FLAG.Plot.sim_step
            plot(sim_data.step.t_soln, sim_data.step.del_phi(:,cell_of_interest),'Linewidth',2,'DisplayName','Step')
        end
        lgn = legend;
        lgn.Location = 'southeast';
        xlabel('Time [s]')
        ylabel('Voltage [V]')
        xlim([xmin,limit])

    end
    if FLAG.Plot.C_Liion
        figure
        hold on
        title('C_{Li^+}')
        if FLAG.Plot.Slink.plant_noise
            plot(t_slink_plant_noise, y_slink_plant_noise(:,localP.C_Liion),'Linewidth',2,'DisplayName','Slink Plant Noise')
        end
        if FLAG.Plot.Slink.ROM
            plot(t_slink_ROM        , y_slink_ROM(:,localP.C_Liion),'o','Linewidth',2,'DisplayName','Slink ROM')
        end
        if FLAG.Plot.Slink.Est
            plot(t_slink_est(start_idx_est:end), y_slink_est(start_idx_est:end,localP.C_Liion),'o','Linewidth',2,'DisplayName','Slink Est')
        end
        if FLAG.Plot.Slink.Est_asy
            plot(t_slink_est_asy(start_idx_est_asy:end), y_slink_est_asy(start_idx_est_asy:end,localP.C_Liion),'o','Linewidth',2,'DisplayName','Slink Est Asy')
        end        
        if FLAG.Plot.Slink.plant
            plot(t_slink_plant      , y_slink_plant(:,localP.C_Liion),'Linewidth',2,'DisplayName','Slink Plant')
        end
        if FLAG.Plot.sim_impl
            plot(sim_data.impl.t_soln, sim_data.impl.C_Liion(:,cell_of_interest),'Linewidth',2,'DisplayName','Impluse')
        end
        if FLAG.Plot.sim_sine
            plot(sim_data.sine.t_soln, sim_data.sine.C_Liion(:,cell_of_interest),'Linewidth',2,'DisplayName','Harmonic')
        end
        if FLAG.Plot.ROM_step
            plot(t_SS_ROM,y_SS_ROM(:,localP.C_Liion)+IC(localP.C_Liion),'o','Linewidth',2,'DisplayName','ROM Step')
        end
        if FLAG.Plot.sim_step
            plot(sim_data.step.t_soln, sim_data.step.C_Liion(:,cell_of_interest),'Linewidth',2,'DisplayName','Step')
        end
        lgn = legend;
        lgn.Location = 'southeast';
        xlabel('Time [s]')
        ylabel('Concentration [kmol m^{-3}]')
        xlim([xmin,limit])

    end
    if FLAG.Plot.X_surf
        figure
        hold on
        title('X_{surf}')
        if FLAG.Plot.Slink.plant_noise
            plot(t_slink_plant_noise, y_slink_plant_noise(:,localP.X_Li_surf),'Linewidth',2,'DisplayName','Slink Plant Noise')
        end
        if FLAG.Plot.Slink.ROM
            plot(t_slink_ROM        , y_slink_ROM(:,localP.X_Li_surf),'o','Linewidth',2,'DisplayName','Slink ROM')
        end
        if FLAG.Plot.Slink.Est
            plot(t_slink_est(start_idx_est:end), y_slink_est(start_idx_est:end,localP.X_Li_surf),'o','Linewidth',2,'DisplayName','Slink Est')
        end
        if FLAG.Plot.Slink.Est_asy
            plot(t_slink_est_asy(start_idx_est_asy:end), y_slink_est_asy(start_idx_est_asy:end,localP.X_Li_surf),'o','Linewidth',2,'DisplayName','Slink Est Asy')
        end        
        if FLAG.Plot.Slink.plant
            plot(t_slink_plant      , y_slink_plant(:,localP.X_Li_surf),'Linewidth',2,'DisplayName','Slink Plant')
        end
        if FLAG.Plot.sim_impl
            plot(sim_data.impl.t_soln, sim_data.impl.X_Li_surf(:,cell_of_interest),'Linewidth',2,'DisplayName','Impluse')
        end
        if FLAG.Plot.sim_sine
            plot(sim_data.sine.t_soln, sim_data.sine.X_Li_surf(:,cell_of_interest),'Linewidth',2,'DisplayName','Harmonic')
        end
        if FLAG.Plot.ROM_step
            plot(t_SS_ROM,y_SS_ROM(:,localP.X_Li_surf)+IC(localP.X_Li_surf),'o','Linewidth',2,'DisplayName','ROM Step')
        end
        if FLAG.Plot.sim_step
            plot(sim_data.step.t_soln, sim_data.step.X_Li_surf(:,cell_of_interest),'Linewidth',2,'DisplayName','Step')
        end
        lgn = legend;
        xlabel('Time [s]')
        ylabel('Mole Fraction [-]')
        xlim([xmin,limit])
        
    end
    if FLAG.Plot.i_Far
        figure
        hold on
        title('i_{Far}')
        if FLAG.Plot.Slink.plant_noise
            plot(t_slink_plant_noise, y_slink_plant_noise(:,localP.i_Far),'Linewidth',2,'DisplayName','Slink Plant Noise')
        end
        if FLAG.Plot.Slink.ROM
            plot(t_slink_ROM        , y_slink_ROM(:,localP.i_Far),'o','Linewidth',2,'DisplayName','Slink ROM')
        end
        if FLAG.Plot.Slink.Est
            plot(t_slink_est(start_idx_est:end), y_slink_est(start_idx_est:end,localP.i_Far),'o','Linewidth',2,'DisplayName','Slink Est')
        end
        if FLAG.Plot.Slink.Est_asy
            plot(t_slink_est_asy(start_idx_est_asy:end), y_slink_est_asy(start_idx_est_asy:end,localP.i_Far),'o','Linewidth',2,'DisplayName','Slink Est Asy')
        end        
        if FLAG.Plot.Slink.plant
            plot(t_slink_plant      , y_slink_plant(:,localP.i_Far),'Linewidth',2,'DisplayName','Slink Plant')
        end
        if FLAG.Plot.sim_impl
            plot(sim_data.impl.t_soln, sim_data.impl.i_Far(:,cell_of_interest),'Linewidth',2,'DisplayName','Impluse')
        end
        if FLAG.Plot.sim_sine
            plot(sim_data.sine.t_soln, sim_data.sine.i_Far(:,cell_of_interest),'Linewidth',2,'DisplayName','Harmonic')
        end
        if FLAG.Plot.ROM_step
            plot(t_SS_ROM,y_SS_ROM(:,localP.i_Far)+IC(localP.i_Far),'o','Linewidth',2,'DisplayName','ROM Step')
        end
        if FLAG.Plot.sim_step
            plot(sim_data.step.t_soln, sim_data.step.i_Far(:,cell_of_interest),'Linewidth',2,'DisplayName','Step')
        end
        lgn = legend;
        lgn.Location = 'southeast';
        xlabel('Time [s]')
        ylabel('Current Density [A m^{-2}]')
        xlim([xmin,limit])
        
    end
    if FLAG.Plot.eta
        figure
        hold on
        title('\eta')
        if FLAG.Plot.Slink.plant_noise
            plot(t_slink_plant_noise, y_slink_plant_noise(:,localP.eta),'Linewidth',2,'DisplayName','Slink Plant Noise')
        end
        if FLAG.Plot.Slink.ROM
            plot(t_slink_ROM        , y_slink_ROM(:,localP.eta),'o','Linewidth',2,'DisplayName','Slink ROM')
        end
        if FLAG.Plot.Slink.Est
            plot(t_slink_est(start_idx_est:end), y_slink_est(start_idx_est:end,localP.eta),'o','Linewidth',2,'DisplayName','Slink Est')
        end
        if FLAG.Plot.Slink.Est_asy
            plot(t_slink_est_asy(start_idx_est_asy:end), y_slink_est_asy(start_idx_est_asy:end,localP.eta),'o','Linewidth',2,'DisplayName','Slink Est Asy')
        end        
        if FLAG.Plot.Slink.plant
            plot(t_slink_plant      , y_slink_plant(:,localP.eta),'Linewidth',2,'DisplayName','Slink Plant')
        end
        if FLAG.Plot.sim_impl
            plot(sim_data.impl.t_soln, sim_data.impl.eta(:,cell_of_interest),'Linewidth',2,'DisplayName','Impluse')
        end
        if FLAG.Plot.sim_sine
            plot(sim_data.sine.t_soln, sim_data.sine.eta(:,cell_of_interest),'Linewidth',2,'DisplayName','Harmonic')
        end
        if FLAG.Plot.ROM_step
            plot(t_SS_ROM,y_SS_ROM(:,localP.eta)+IC(localP.eta),'o','Linewidth',2,'DisplayName','ROM Step')
        end
        if FLAG.Plot.sim_step
            plot(sim_data.step.t_soln, sim_data.step.eta(:,cell_of_interest),'Linewidth',2,'DisplayName','Step')
        end
        lgn = legend;
        lgn.Location = 'southeast';
        xlabel('Time [s]')
        ylabel('Voltage [V]')
        xlim([xmin,limit])
        
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

end