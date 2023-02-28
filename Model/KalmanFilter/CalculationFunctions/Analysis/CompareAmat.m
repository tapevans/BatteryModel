%% Compare A matrix
% Want to compare the eigenvalues of A generated from the ROM across a
% range of SOC
clear all; close all; clc;

%% Parameters
SOC_center = 50; % Centered SOC
SOC_pm = 7; % Number plus/minus around it
Ts = 10; % Sampling Time
N_compare = 10; % How many eigen values to compare


%% Create SOC vec
SOC_min = SOC_center - SOC_pm;
SOC_max = SOC_center + SOC_pm;
if SOC_min < 0
    SOC_min = 0;
end
if SOC_max > 100
    SOC_max = 100;
end

SOC_vec = SOC_min:1:SOC_max;
N_SOC = length(SOC_vec);

FLAG.OverwriteData = 0;

%% Load Results
    RUNSIM = false;
    save_filename = ['SOC_min' num2str(SOC_min) '_SOCmax' num2str(SOC_max) '_Ts' num2str(Ts) '.mat'];
    if isfile(save_filename)
        if FLAG.OverwriteData
            disp('Overwriting Data File')
            RUNSIM = true;
            delete(save_filename)
        else
            disp('File Exists')
            load(save_filename);
        end
    else
        RUNSIM = true;
    end

if RUNSIM
    %% N,P,RESULTS
    % Can all be empty since I'm not plotting impluse response
    N = struct();
    P = struct();
    RESULTS = struct();


    %% FLAG
    FLAG.Ts = Ts;
    FLAG.Analysis.PlotImp = 0;
    FLAG.folderpath = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\ObservabilityTest';


    %% Desired Outputs
    FLAG.DesOut.cell_voltage  = 1; % Terminal Voltage
    FLAG.DesOut.delta_phi     = 1; % Electrostatic potential difference between active material and electrolyte @AN/SEP interface
    FLAG.DesOut.i_Far         = 1; % Chemical reaction current density at SEI @AN/SEP interface
    FLAG.DesOut.eta           = 1; % Overpotential at SEI @AN/SEP interface
    FLAG.DesOut.C_Liion       = 1; % Concentration of Li^+ in the electrolyte @AN/SEP interface
    FLAG.DesOut.C_Li          = 1; % Concentration of Li at the surface of active material @AN/SEP interface
    FLAG.DesOut.delta_C_Li    = 1; % Difference of Li concentration between the surface and center of active material particle @AN/SEP interface
    FLAG.DesOut.T             = 1; % Temperature of control volume @AN/SEP interface


    %% Set pointers for variables
    i = 1;
    if FLAG.DesOut.cell_voltage   % Cell Voltage
        P.cell_voltage = i;
        RESULTS.Labels.title{i} = 'Cell Voltage';
        RESULTS.Labels.unit{i}  = 'Voltage [V]';
        i = i + 1;
    end
    if FLAG.DesOut.delta_phi      % Delta Phi
        P.delta_phi = i;
        RESULTS.Labels.title{i} = '\Delta \phi';
        RESULTS.Labels.unit{i}  = 'Voltage [V]';
        i = i + 1;
    end
    if FLAG.DesOut.i_Far          % i_Far
        P.i_Far = i;
        RESULTS.Labels.title{i} = 'i_{Far}';
        RESULTS.Labels.unit{i}  = 'Current [A/m^2]';
        i = i + 1;
    end
    if FLAG.DesOut.eta          % eta
        P.eta = i;
        RESULTS.Labels.title{i} = '\eta';
        RESULTS.Labels.unit{i}  = 'Voltage [V]';
        i = i + 1;
    end
    if FLAG.DesOut.C_Liion        % C_Li+
        P.C_Liion = i;
        RESULTS.Labels.title{i} = 'C_{Li^+}';
        RESULTS.Labels.unit{i}  = 'Concentration [kmol/m^3]';
        i = i + 1;
    end
    if FLAG.DesOut.C_Li           % C_Li
        P.C_Li = i;
        RESULTS.Labels.title{i} = 'C_{Li,surf}';
        RESULTS.Labels.unit{i}  = 'Concentration [kmol/m^3]';
        i = i + 1;
    end
    if FLAG.DesOut.delta_C_Li     % Delta C_Li
        P.delta_C_Li = i;
        RESULTS.Labels.title{i} = '\Delta C_{Li}';
        RESULTS.Labels.unit{i}  = 'Concentration [kmol/m^3]';
        i = i + 1;
    end
    if FLAG.DesOut.T              % Temperature
        P.T = i;
        RESULTS.Labels.title{i} = 'Temperature';
        RESULTS.Labels.unit{i}  = 'Temperature [K]';
        i = i + 1;
    end

    N.DesOut = i-1;


    %% Output Matrix
    % Load a file to get pointers
    FLAG.SOC = SOC_center;
    filename = getImpulseFilename(FLAG);
    %     filename = getSSFilename(FLAG);
    sys = load(filename);

    N.N_In  = 1;
    % I_user
    % Outputs
    % Cell Voltage
    % Delta Phi   @AN/SEP
    % Temperature @AN/SEP
    % C_Liion     @AN/SEP
    % X_surf      @AN/SEP
    % i_Far       @AN/SEP
    % eta         @AN/SEP
    SIM.OutputMatrix = zeros(N.DesOut , sys.N.N_SV_tot);
    if FLAG.DesOut.cell_voltage
        % Cell Voltage
        idx_phi_ed_AN = sys.P.phi_ed;

        i = sys.N.N_CV_CA(end);
        index_offset = (i-1)*sys.N.N_SV_CA + sys.N.N_SV_AN_tot + sys.N.N_SV_SEP_tot;
        idx_phi_ed_CA = index_offset + sys.P.phi_ed;

        SIM.OutputMatrix(P.cell_voltage,idx_phi_ed_AN) = -1;
        SIM.OutputMatrix(P.cell_voltage,idx_phi_ed_CA) =  1;
    end
    % @AN/SEP
    i = sys.N.N_CV_AN(end);
    index_offset = (i-1)*sys.N.N_SV_AN;
    if FLAG.DesOut.delta_phi
        % Delta Phi   @AN/SEP
        idx = index_offset + sys.P.del_phi;
        SIM.OutputMatrix(P.delta_phi,idx) =  1;
    end
    if FLAG.DesOut.i_Far
        % i_Far      @AN/SEP
        idx = index_offset + sys.P.i_PS;
        SIM.OutputMatrix(P.i_Far,idx) = 1;
    end
    if FLAG.DesOut.eta
        % eta       @AN/SEP
        idx = index_offset + sys.P.V_2;
        SIM.OutputMatrix(P.eta,idx) = 1;
        idx = index_offset + sys.P.V_1;
        SIM.OutputMatrix(P.eta,idx) = -1;
    end
    if FLAG.DesOut.C_Liion
        % C_Liion     @AN/SEP
        idx = index_offset + sys.P.C_Liion;
        SIM.OutputMatrix(P.C_Liion,idx) = 1;
    end
    if FLAG.DesOut.C_Li
        % C_Li,surf   @AN/SEP
        idx = index_offset + sys.P.C_Li_surf_AN;
        SIM.OutputMatrix(P.C_Li,idx) = 1;%/AN.C_Li_max;
    end
    if FLAG.DesOut.delta_C_Li
        % Delta C_Li  @AN/SEP (Over entire radius of particle)
        idx = index_offset + sys.P.C_Li_surf_AN; % Surface Node
        SIM.OutputMatrix(P.delta_C_Li,idx) = 1;
        idx = index_offset + sys.P.C_Li;         % Most Interior Node
        SIM.OutputMatrix(P.delta_C_Li,idx) = -1;
    end
    if FLAG.DesOut.T
        % Temperature @AN/SEP
        idx = index_offset + sys.P.T;
        SIM.OutputMatrix(P.T,idx) = 1;
    end


    %% Loop through all SOC
    for SS = 1:N_SOC
        FLAG.SOC = SOC_vec(SS);
        SIM.Ts = Ts;
        [ROM_sys{SS}] = getHoKalmanROM(SIM,N,P,FLAG,RESULTS);
    end


    %% Loop through all ROM and save the first N_compare eigienvalues of A
    eigen_mat = nan(N_SOC,N_compare);
    for SS = 1:N_SOC
        eigen_val = abs( real( eig( ROM_sys{SS}.A) ) );
        eigen_val = sort(eigen_val,'descend');
        eigen_mat(SS,:) = eigen_val(1:N_compare)';
    end

    %% Save Results
    save(save_filename,'eigen_mat','SOC_min','SOC_max','SOC_vec','Ts')

end


%% Plot Results
    [~,N_lambda] = size(eigen_mat);
    [X,Y] = meshgrid(1:1:N_lambda  , SOC_vec);
    f = figure;
    [~,h] = contourf(X,Y,eigen_mat);
    title(['Comparing the first ' num2str(N_lambda) ' of the ROM at T_s=' num2str(Ts)])
    xlabel('\lambda_i')
    ylabel('SOC [%]')
    colorbar
%     colormap(f, flipud(colormap(colormapping)))
    clim([0.5 1])
    h.LineStyle = 'none';
    pic_filename = save_filename(1:end-4);
    exportgraphics(gca,[pic_filename '.png'],'Resolution',500)









