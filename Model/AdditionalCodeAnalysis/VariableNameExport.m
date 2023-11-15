%% Output SV Names to CSV with Output Matrix
clear all;
close all; 
clc;

%%
filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestImpedanceContributions\Standard_SS_EIS_SOC50.mat';
load(filename)


%% Pointers
i = 1;
P.Name      = i; i = i+1;
P.CV_Number = i; i = i+1;
P.SV_Number = i; i = i+1;
P.ANCV      = i; i = i+1;
P.SEPCV     = i; i = i+1;
P.CACV      = i; i = i+1;

NumP = i-1;

%%
Data = cell(N.N_SV_tot , NumP);

% ---- Anode ----
    for i = 1:N.N_CV_AN
        index_offset = (i-1)*N.N_SV_AN;

        CV_Number = i;
        ANCV      = i;
        SEPCV     = '';
        CACV      = '';

        % Temp
            Data{index_offset+P.T       , P.Name }      =  'Temperature'; 
            Data{index_offset+P.T       , P.CV_Number } =  CV_Number;
            Data{index_offset+P.T       , P.SV_Number } =  index_offset+P.T;
            Data{index_offset+P.T       , P.ANCV }      =  ANCV;
            Data{index_offset+P.T       , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.T       , P.CACV }      =  CACV;
        % del_phi
            Data{index_offset+P.del_phi , P.Name }      =  'del_phi';
            Data{index_offset+P.del_phi , P.CV_Number } =  CV_Number;
            Data{index_offset+P.del_phi , P.SV_Number } =  index_offset+P.del_phi;
            Data{index_offset+P.del_phi , P.ANCV }      =  ANCV;
            Data{index_offset+P.del_phi , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.del_phi , P.CACV }      =  CACV;
        % phi_ed 
            Data{index_offset+P.phi_ed  , P.Name }      =  'phi_ed';
            Data{index_offset+P.phi_ed  , P.CV_Number } =  CV_Number;
            Data{index_offset+P.phi_ed  , P.SV_Number } =  index_offset+P.phi_ed;
            Data{index_offset+P.phi_ed  , P.ANCV }      =  ANCV;
            Data{index_offset+P.phi_ed  , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.phi_ed  , P.CACV }      =  CACV;
        % V_1
            Data{index_offset+P.V_1     , P.Name }      =  'V_1';
            Data{index_offset+P.V_1     , P.CV_Number } =  CV_Number;
            Data{index_offset+P.V_1     , P.SV_Number } =  index_offset+P.V_1;
            Data{index_offset+P.V_1     , P.ANCV }      =  ANCV;
            Data{index_offset+P.V_1     , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.V_1     , P.CACV }      =  CACV;
        % V_2
            Data{index_offset+P.V_2     , P.Name }      =  'V_2';
            Data{index_offset+P.V_2     , P.CV_Number } =  CV_Number;
            Data{index_offset+P.V_2     , P.SV_Number } =  index_offset+P.V_2;
            Data{index_offset+P.V_2     , P.ANCV }      =  ANCV;
            Data{index_offset+P.V_2     , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.V_2     , P.CACV }      =  CACV;
        % i_PS
            Data{index_offset+P.i_PS    , P.Name }      =  'i_PS';
            Data{index_offset+P.i_PS    , P.CV_Number } =  CV_Number;
            Data{index_offset+P.i_PS    , P.SV_Number } =  index_offset+P.i_PS;
            Data{index_offset+P.i_PS    , P.ANCV }      =  ANCV;
            Data{index_offset+P.i_PS    , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.i_PS    , P.CACV }      =  CACV;
        % C_Li^+
            Data{index_offset+P.C_Liion , P.Name }      =  'C_Liion';
            Data{index_offset+P.C_Liion , P.CV_Number } =  CV_Number;
            Data{index_offset+P.C_Liion , P.SV_Number } =  index_offset+P.C_Liion;
            Data{index_offset+P.C_Liion , P.ANCV }      =  ANCV;
            Data{index_offset+P.C_Liion , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.C_Liion , P.CACV }      =  CACV;
        % C_Li
            for j = 1:N.N_R_AN
                Data{index_offset+P.C_Li+j-1 , P.Name}       = ['C_Li_R' num2str(j)];
                Data{index_offset+P.C_Li+j-1 , P.CV_Number } =  CV_Number;
                Data{index_offset+P.C_Li+j-1 , P.SV_Number } =  index_offset+P.C_Li+j-1;
                Data{index_offset+P.C_Li+j-1 , P.ANCV }      =  ANCV;
                Data{index_offset+P.C_Li+j-1 , P.SEPCV }     =  SEPCV;
                Data{index_offset+P.C_Li+j-1 , P.CACV }      =  CACV;
            end
    end

% ---- Separator ----
    for i = 1:N.N_CV_SEP
        index_offset = (i-1)*N.N_SV_SEP + N.N_SV_AN_tot;

        CV_Number = i+N.N_CV_AN;
        ANCV      = '';
        SEPCV     = i;
        CACV      = '';

        % Temp
            Data{index_offset+P.SEP.T       , P.Name }      =  'Temperature';
            Data{index_offset+P.SEP.T       , P.CV_Number } =  CV_Number;
            Data{index_offset+P.SEP.T       , P.SV_Number } =  index_offset+P.SEP.T;
            Data{index_offset+P.SEP.T       , P.ANCV }      =  ANCV;
            Data{index_offset+P.SEP.T       , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.SEP.T       , P.CACV }      =  CACV;
        % phi_el
            Data{index_offset+P.SEP.phi_el  , P.Name }      =  'phi_el';
            Data{index_offset+P.SEP.phi_el  , P.CV_Number } =  CV_Number;
            Data{index_offset+P.SEP.phi_el  , P.SV_Number } =  index_offset+P.SEP.phi_el;
            Data{index_offset+P.SEP.phi_el  , P.ANCV }      =  ANCV;
            Data{index_offset+P.SEP.phi_el  , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.SEP.phi_el  , P.CACV }      =  CACV;
        % C_Li^+
            Data{index_offset+P.SEP.C_Liion , P.Name }      =  'C_Liion'; 
            Data{index_offset+P.SEP.C_Liion , P.CV_Number } =  CV_Number;
            Data{index_offset+P.SEP.C_Liion , P.SV_Number } =  index_offset+P.SEP.C_Liion;
            Data{index_offset+P.SEP.C_Liion , P.ANCV }      =  ANCV;
            Data{index_offset+P.SEP.C_Liion , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.SEP.C_Liion , P.CACV }      =  CACV;
    end

% ---- Cathode ----
    for i = 1:N.N_CV_CA
        index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;

        CV_Number = i + N.N_CV_AN + N.N_CV_SEP;
        ANCV      = '';
        SEPCV     = '';
        CACV      = i;

        % Temp
            Data{index_offset+P.T       , P.Name }      =  'Temperature'; 
            Data{index_offset+P.T       , P.CV_Number } =  CV_Number;
            Data{index_offset+P.T       , P.SV_Number } =  index_offset+P.T;
            Data{index_offset+P.T       , P.ANCV }      =  ANCV;
            Data{index_offset+P.T       , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.T       , P.CACV }      =  CACV;
        % del_phi
            Data{index_offset+P.del_phi , P.Name }      =  'del_phi';
            Data{index_offset+P.del_phi , P.CV_Number } =  CV_Number;
            Data{index_offset+P.del_phi , P.SV_Number } =  index_offset+P.del_phi;
            Data{index_offset+P.del_phi , P.ANCV }      =  ANCV;
            Data{index_offset+P.del_phi , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.del_phi , P.CACV }      =  CACV;
        % phi_ed 
            Data{index_offset+P.phi_ed  , P.Name }      =  'phi_ed';
            Data{index_offset+P.phi_ed  , P.CV_Number } =  CV_Number;
            Data{index_offset+P.phi_ed  , P.SV_Number } =  index_offset+P.phi_ed;
            Data{index_offset+P.phi_ed  , P.ANCV }      =  ANCV;
            Data{index_offset+P.phi_ed  , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.phi_ed  , P.CACV }      =  CACV;
        % V_1
            Data{index_offset+P.V_1     , P.Name }      =  'V_1';
            Data{index_offset+P.V_1     , P.CV_Number } =  CV_Number;
            Data{index_offset+P.V_1     , P.SV_Number } =  index_offset+P.V_1;
            Data{index_offset+P.V_1     , P.ANCV }      =  ANCV;
            Data{index_offset+P.V_1     , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.V_1     , P.CACV }      =  CACV;
        % V_2
            Data{index_offset+P.V_2     , P.Name }      =  'V_2';
            Data{index_offset+P.V_2     , P.CV_Number } =  CV_Number;
            Data{index_offset+P.V_2     , P.SV_Number } =  index_offset+P.V_2;
            Data{index_offset+P.V_2     , P.ANCV }      =  ANCV;
            Data{index_offset+P.V_2     , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.V_2     , P.CACV }      =  CACV;
        % i_PS
            Data{index_offset+P.i_PS    , P.Name }      =  'i_PS';
            Data{index_offset+P.i_PS    , P.CV_Number } =  CV_Number;
            Data{index_offset+P.i_PS    , P.SV_Number } =  index_offset+P.i_PS;
            Data{index_offset+P.i_PS    , P.ANCV }      =  ANCV;
            Data{index_offset+P.i_PS    , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.i_PS    , P.CACV }      =  CACV;
        % C_Li^+
            Data{index_offset+P.C_Liion , P.Name }      =  'C_Liion';
            Data{index_offset+P.C_Liion , P.CV_Number } =  CV_Number;
            Data{index_offset+P.C_Liion , P.SV_Number } =  index_offset+P.C_Liion;
            Data{index_offset+P.C_Liion , P.ANCV }      =  ANCV;
            Data{index_offset+P.C_Liion , P.SEPCV }     =  SEPCV;
            Data{index_offset+P.C_Liion , P.CACV }      =  CACV;
        % C_Li
            for j = 1:N.N_R_CA
                Data{index_offset+P.C_Li+j-1 , P.Name}       = ['C_Li_R' num2str(j)];
                Data{index_offset+P.C_Li+j-1 , P.CV_Number } =  CV_Number;
                Data{index_offset+P.C_Li+j-1 , P.SV_Number } =  index_offset+P.C_Li+j-1;
                Data{index_offset+P.C_Li+j-1 , P.ANCV }      =  ANCV;
                Data{index_offset+P.C_Li+j-1 , P.SEPCV }     =  SEPCV;
                Data{index_offset+P.C_Li+j-1 , P.CACV }      =  CACV;
            end
    end