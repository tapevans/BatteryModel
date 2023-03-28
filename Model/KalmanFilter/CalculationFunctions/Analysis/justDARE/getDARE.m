function [RESULTS] = getDARE(SIM,FLAG,N,P,RESULTS)
%% Initialize
    % Qmin = -7;
    % Qmax =  2;
    % Qnum = 20;
    Qmin = -7;
    Qmax =  5;
    Qnum =  13;%20
    Q_vec = logspace(Qmin , Qmax , Qnum);
    %Q_vec = [1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2];

    % Rmin = -7;
    % Rmax =  2;
    % Rnum = 20;
    Rmin = -7;
    Rmax =  4;
    Rnum =  23;%20
    R_vec = logspace(Rmin , Rmax , Rnum);
    %R_vec = [1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2];
    
    RESULTS.DARE.Q_vec = Q_vec;
    RESULTS.DARE.R_vec = R_vec;
    RESULTS.DARE.IDV.CPCT = nan(Qnum , Rnum , N.DesOut);
    RESULTS.DARE.COM.CPCT = nan(Qnum , Rnum , N.DesOut);


%% Get ROM
    %  1) Matlab SS_DT
    %  2) Ho-Kalman
    switch FLAG.EstimatorModel
        case 1
            % [~ , sys_DT] = getSS_System(SIM,N,P,FLAG);
            % est_sys_tot{1} = sys_DT;
        case 2
            [HK_sys] = getHoKalmanROM(SIM,N,P,FLAG,RESULTS);
            est_sys_tot = HK_sys;
    end


%% Loop through all Q, R, and outputs
% Individual
if FLAG.IDV
    for RR = 1:Rnum
        % RR
        SIM.R_0 = R_vec(RR);
        for QQ = 1:Qnum
            % QQ
            SIM.Qi = Q_vec(QQ);
            for OO = 1:N.DesOut
                if RR == 17 && QQ == 5 && OO == 4
                    R = R_vec(RR) 
                    Q = Q_vec(QQ)
                    RESULTS.Labels.title{OO}
                end
                est_sys = est_sys_tot{OO};
                [~, P_inf] = AsymptoticPreCalcs(FLAG,SIM,est_sys);
                CPCT = est_sys.C * P_inf * est_sys.C';
                RESULTS.DARE.IDV.CPCT(QQ , RR , OO) = CPCT(2,2);
            end
        end
    end
end

% Combined
if FLAG.COM
    for RR = 1:Rnum
        SIM.R_0 = R_vec(RR);
        for QQ = 1:Qnum
            SIM.Qi = Q_vec(QQ);

            est_sys = est_sys_tot{end};
            [~, P_inf] = AsymptoticPreCalcs(FLAG,SIM,est_sys);
            CPCT = est_sys.C * P_inf * est_sys.C';
            RESULTS.DARE.COM.CPCT(QQ , RR , :) = diag(CPCT);
        end
    end
end

%% Plot Data
if FLAG.IDV
    Z = log10(RESULTS.DARE.IDV.CPCT);
    Q_mag_vec = linspace(Qmin , Qmax , Qnum);
    R_mag_vec = linspace(Rmin , Rmax , Rnum);
    [X,Y] = meshgrid(R_mag_vec , Q_mag_vec);
    for OO = 1:N.DesOut
        figure
        s = surf(X,Y,Z(:,:,OO));
        s.LineStyle = 'none';
        xlabel('R')
        ylabel('Q')
        zlabel('Error Covar (CPC^T)')
        title([RESULTS.Labels.title{OO} ' Error Covariance Individual'])
        xlim([Rmin,Rmax])
        ylim([Qmin,Qmax])
        % clim([0,200])
        % colorbar
    end
end

if FLAG.COM
    Z = log10(RESULTS.DARE.COM.CPCT);
    Q_mag_vec = linspace(Qmin , Qmax , Qnum);
    R_mag_vec = linspace(Rmin , Rmax , Rnum);
    [X,Y] = meshgrid(R_mag_vec , Q_mag_vec);
    for OO = 1:N.DesOut
        figure
        s = surf(X,Y,Z(:,:,OO));
        s.LineStyle = 'none';
        xlabel('R')
        ylabel('Q')
        zlabel('Error Covar (CPC^T)')
        title([RESULTS.Labels.title{OO} ' Error Covariance Combined'])
        xlim([Rmin,Rmax])
        ylim([Qmin,Qmax])
        % clim([0,200])
        % colorbar
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