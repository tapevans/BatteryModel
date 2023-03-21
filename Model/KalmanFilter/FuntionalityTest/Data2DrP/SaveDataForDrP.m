%% Save Data for Dr. P
clear all; close all; clc;

%% Parameters
FLAG.PRBSAmp  = 1;
FLAG.SOC      = 50;
FLAG.Tswitch  = 100;
FLAG.Ts       = 1;
FLAG.UseRelax = 1;

%% Load Data
if FLAG.UseRelax
    folderpathPRBS     = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\zeroMeanPRBS_withRelax_Sims';
else
    folderpathPRBS     = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\PRBS_Sims';
end
just_file = ['PRBS_Sims_PRBS_Amp' num2str(FLAG.PRBSAmp) '_SOC' num2str(FLAG.SOC) '_SwitchingTime' num2str(FLAG.Tswitch) '.mat'];
simData = load([folderpathPRBS filesep just_file]);


%% Get Pointers
i = 1;
    %if FLAG.DesOut.cell_voltage   % Cell Voltage
        P.cell_voltage = i; 
        Labels.title{i} = 'Cell Voltage';
        Labels.unit{i}  = 'Voltage [V]';
        %foldername{i}   = 'CellVoltage';
        i = i + 1;       
    %end
    %if FLAG.DesOut.delta_phi      % Delta Phi
        P.delta_phi = i; 
        Labels.title{i} = '\Delta \phi';
        Labels.unit{i}  = 'Voltage [V]';
        %foldername{i}   = 'DeltaPhi';
        i = i + 1;    
    %end
%     if FLAG.DesOut.i_Far          % i_Far
%         P.i_Far = i;
%         Labels.title{i} = 'i_{Far}';
%         Labels.unit{i}  = 'Current [A/m^2]';
%         foldername{i}   = 'iFar';
%         i = i + 1;
%     end
    %if FLAG.DesOut.eta          % eta
        P.eta = i;
        Labels.title{i} = '\eta';
        Labels.unit{i}  = 'Voltage [V]';
        %foldername{i}   = 'eta';
        i = i + 1;
    %end
    %if FLAG.DesOut.C_Liion        % C_Li+
        P.C_Liion = i; 
        Labels.title{i} = 'C_{Li^+}';
        Labels.unit{i}  = 'Concentration [kmol/m^3]';
        %foldername{i}   = 'C_Liion';
        i = i + 1;
    %end
    %if FLAG.DesOut.C_Li           % C_Li
        P.C_Li = i;
        Labels.title{i} = 'C_{Li,surf}';
        Labels.unit{i}  = 'Concentration [kmol/m^3]';
        %foldername{i}   = 'C_Li_surf';
        i = i + 1;
    %end
    %if FLAG.DesOut.delta_C_Li     % Delta C_Li
        P.delta_C_Li = i;
        Labels.title{i} = '\Delta C_{Li}';
        Labels.unit{i}  = 'Concentration [kmol/m^3]';
        %foldername{i}   = 'DeltaC_Li';
        i = i + 1;
    %end
%     if FLAG.DesOut.T              % Temperature
%         P.T = i;
%         Labels.title{i} = 'Temperature';
%         Labels.unit{i}  = 'Temperature [K]';
%         foldername{i}   = 'Temperature';
%         i = i + 1;
%     end

    N.DesOut = i-1;

%% Output Matrix
    N.N_In  = 1;
        % i_user

    % Outputs
        % Cell Voltage
        % Delta Phi   @AN/SEP
        % i_Far       @AN/SEP
        % eta         @AN/SEP
        % C_Liion     @AN/SEP
        % C_Li        @AN/SEP
        % Delta C_Li  @AN/SEP
        % Temperature @AN/SEP
    SIM.OutputMatrix = zeros(N.DesOut , simData.N.N_SV_tot);
        %if FLAG.DesOut.cell_voltage
        % Cell Voltage
            idx_phi_ed_AN = simData.P.phi_ed;

            i = simData.N.N_CV_CA(end);
            index_offset = (i-1)*simData.N.N_SV_CA + simData.N.N_SV_AN_tot + simData.N.N_SV_SEP_tot;
            idx_phi_ed_CA = index_offset + simData.P.phi_ed;

            SIM.OutputMatrix(P.cell_voltage,idx_phi_ed_AN) = -1; 
            SIM.OutputMatrix(P.cell_voltage,idx_phi_ed_CA) =  1;
        %end
        % @AN/SEP
            i = simData.N.N_CV_AN(end);
            index_offset = (i-1)*simData.N.N_SV_AN;
        %if FLAG.DesOut.delta_phi
        % Delta Phi   @AN/SEP
            idx = index_offset + simData.P.del_phi;
            SIM.OutputMatrix(P.delta_phi,idx) =  1; 
        %end
%         %if FLAG.DesOut.i_Far
%         % i_Far      @AN/SEP
%             idx = index_offset + simData.P.i_PS;
%             SIM.OutputMatrix(P.i_Far,idx) = 1;  
%         %end
        %if FLAG.DesOut.eta
        % eta       @AN/SEP
            idx = index_offset + simData.P.V_2;
            SIM.OutputMatrix(P.eta,idx) = 1;
            idx = index_offset + simData.P.V_1;
            SIM.OutputMatrix(P.eta,idx) = -1;
        %end
        %if FLAG.DesOut.C_Liion
        % C_Liion     @AN/SEP
            idx = index_offset + simData.P.C_Liion;
            SIM.OutputMatrix(P.C_Liion,idx) = 1;
        %end
        %if FLAG.DesOut.C_Li
        % C_Li,surf   @AN/SEP
            idx = index_offset + simData.P.C_Li_surf_AN;
            SIM.OutputMatrix(P.C_Li,idx) = 1;%/AN.C_Li_max;
        %end
        %if FLAG.DesOut.delta_C_Li
        % Delta C_Li  @AN/SEP (Over entire radius of particle)
            idx = index_offset + simData.P.C_Li_surf_AN; % Surface Node
            SIM.OutputMatrix(P.delta_C_Li,idx) = 1;
            idx = index_offset + simData.P.C_Li;         % Most Interior Node
            SIM.OutputMatrix(P.delta_C_Li,idx) = -1;
        %end
%         %if FLAG.DesOut.T
%         % Temperature @AN/SEP
%             idx = index_offset + simData.P.T;
%             SIM.OutputMatrix(P.T,idx) = 1; 
%         %end


%% Save Time Data
clear idx
t_DT_vec = 0:FLAG.Ts:simData.t_soln(end);
for i = 1:length(t_DT_vec)
    [~,idx(i)] = min(abs(t_DT_vec(i)-simData.t_soln));
end

SV_soln = simData.SIM.OutputMatrix * simData.SV_soln';

time = t_DT_vec;
u    = simData.i_user(idx)';

outputs(P.cell_voltage,:) = SV_soln(simData.P.OM.cell_volt,idx);
outputs(P.delta_phi   ,:) = SV_soln(simData.P.OM.delta_phi,idx);
outputs(P.eta         ,:) = SV_soln(simData.P.OM.eta,idx);
outputs(P.C_Liion     ,:) = SV_soln(simData.P.OM.C_Liion,idx);
outputs(P.C_Li        ,:) = SV_soln(simData.P.OM.C_Li,idx);
outputs(P.delta_C_Li  ,:) = SV_soln(simData.P.OM.delta_C_Li,idx);


%% Test Plot
figure
plot(time,u,'o')

figure
plot(time,outputs(P.cell_voltage,:))

%% Save Data
if FLAG.UseRelax
    filename = ['PRBS_Relax_Amp' num2str(FLAG.PRBSAmp) '_SOC' num2str(FLAG.SOC) '_SwitchingTime' num2str(FLAG.Tswitch) '_Ts' num2str(FLAG.Ts) '.mat'];
else
    filename = ['PRBS_Amp' num2str(FLAG.PRBSAmp) '_SOC' num2str(FLAG.SOC) '_SwitchingTime' num2str(FLAG.Tswitch) '_Ts' num2str(FLAG.Ts) '.mat'];
end
save(filename,'time','u','outputs','P','Labels')
