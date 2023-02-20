%% 
function Kalman_plotFnc(RESULTS,N)
%%
FLAG.PLOT.NoNoiseCompare = 1;    
    FLAG.PLOT.ode        = 1;
    FLAG.PLOT.SS_CT      = 0;
    FLAG.PLOT.SS_DT      = 0;
    FLAG.PLOT.ROM_HoKal  = 1;


%% No Noise Plots
if FLAG.PLOT.NoNoiseCompare
    if isfield(RESULTS,'ROM_HoKal') || isfield(RESULTS,'SS_DT') || isfield(RESULTS,'SS_CT') || isfield(RESULTS,'ode')
        for i = 1:N.DesOut
            figure
            hold on
            if FLAG.PLOT.ROM_HoKal && isfield(RESULTS,'ROM_HoKal')
                plot(RESULTS.ROM_HoKal.t_soln,RESULTS.ROM_HoKal.z_soln(:,i),'o','LineWidth',2,'DisplayName','ROM Ho-Kal')
            end
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
            lgn.Location = 'best';
            title([RESULTS.Labels.title{i} ' No Noise Comparison'])
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