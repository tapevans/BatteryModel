%% Magnitude of CPCT vs SOC and Ts
%% Define Noise to compare
Q_0 = 1e-6; R_0 = 1e-6;


%% Local FLAGS
localFLAG.Plot.NotNormalized          = 1;
localFLAG.Plot.Normalized2CellVoltage = 1;
localFLAG.Plot.Normalized2MaxCPCT     = 1;

localFLAG.Plot.cell_voltage  = 1; % Terminal Voltage
localFLAG.Plot.delta_phi     = 1; % Electrostatic potential difference between active material and electrolyte @AN/SEP interface
localFLAG.Plot.i_Far         = 1; % Chemical reaction current density at SEI @AN/SEP interface
localFLAG.Plot.eta           = 1; % Overpotential at SEI @AN/SEP interface
localFLAG.Plot.C_Liion       = 1; % Concentration of Li^+ in the electrolyte @AN/SEP interface
localFLAG.Plot.C_Li          = 1; % Concentration of Li at the surface of active material @AN/SEP interface
localFLAG.Plot.delta_C_Li    = 1; % Difference of Li concentration between the surface and center of active material particle @AN/SEP interface
localFLAG.Plot.T             = 1; % Temperature

localFLAG.Overwrite.MyData = 0;

%% Make string for filenames to read in


%% Check if Data has already been consolidated
if % if yes, load it
    if localFLAG.Overwrite.MyData

    end
else % if no, consolidate it

end


    

%% Read in Data
% State of Charge
    SOC_vec = 0:1:100;

% Desired Sampling Rate
    N_t_s   = 25;   % **Keep this fix for now
    T_s_min = -1; % T_s = 10^(T_s_min), **Keep this fix for now
    T_s_max =  1; % T_s = 10^(T_s_max), **Keep this fix for now
    Ts_vec = logspace(T_s_min,T_s_max,N_t_s);


%% Initialize variables for Contour plots


%% Match MyData idx with plot location


%% Plots
% Contour not normalized



% Contour normalized to Cell Voltage



% Countour normalized to max CPCT value

