%% 
function Kalman_plotFnc(RESULTS)
%%
FLAG.PLOT.NoNoiseCompare = 1;
    output_idx = [1 2 3 4 5 6 7];
%     output_idx = [1 2 4 5 6 7];
    output_name = {'Cell Voltage' , '\Delta \phi' , 'Temperature'     , 'C_{Li^+}'                 , 'C_{Li,surf}'              , 'i_{Far}'        , '\eta'}; % '\Delta C_{Li}' ,
    output_unit = {'Voltage [V]'  , 'Voltage [V]' , 'Temperature [K]' , 'Concentration [kmol/m^3]' , 'Concentration [kmol/m^3]' ,'Current [A/m^2]' , 'Voltage [V]'};% 'Concentration [kmol/m^3]' , 
    FLAG.PLOT.ode        = 1;
    FLAG.PLOT.SS_CT      = 0;
    FLAG.PLOT.SS_DT      = 0;
    %FLAG.PLOT.ROM_Mlab   = 0;
    FLAG.PLOT.ROM_HoKal  = 0;




%% No Noise Plots
if FLAG.PLOT.NoNoiseCompare
    if isfield(RESULTS,'ROM_HoKal') || isfield(RESULTS,'ROM_Mlab') || isfield(RESULTS,'SS_DT') || isfield(RESULTS,'SS_CT') || isfield(RESULTS,'ode')
        for i = output_idx
            figure
            hold on
            if FLAG.PLOT.ROM_HoKal && isfield(RESULTS,'ROM_HoKal')
                plot(RESULTS.ROM_HoKal.t_soln,RESULTS.ROM_HoKal.z_soln(:,i),'o','LineWidth',2,'DisplayName','ROM Ho-Kal')
            end
%             if FLAG.PLOT.ROM_Mlab && isfield(RESULTS,'ROM_Mlab')
%                 plot(RESULTS.ROM_Mlab.t_soln,RESULTS.ROM_Mlab.z_soln(:,i),'o','LineWidth',2,'DisplayName','ROM Mlab')
%             end
            if FLAG.PLOT.SS_DT    && isfield(RESULTS,'SS_DT')
                plot(RESULTS.SS_DT.t_soln,RESULTS.SS_DT.z_soln(:,i),'^','LineWidth',2,'DisplayName','SS DT')
            end
            if FLAG.PLOT.SS_CT    && isfield(RESULTS,'SS_CT')
                plot(RESULTS.SS_CT.t_soln,RESULTS.SS_CT.z_soln(:,i),'square','LineWidth',2,'DisplayName','SS CT')
            end
            if FLAG.PLOT.ode      &&  isfield(RESULTS,'ode')
                plot(RESULTS.ode.t_soln,RESULTS.ode.z_soln(i,:),'k','LineWidth',2,'DisplayName','ODE')
            end
            lgn = legend;
            title([output_name{i} ' No Noise Comparison'])
            xlabel('Time [s]')
            ylabel([output_unit{i}])
    
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