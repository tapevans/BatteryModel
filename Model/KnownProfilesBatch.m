clear all; close all; clc;
% I supressed KBSOC and clear all; close all; clc; in CreateProject
% Then call this to loop through each SOC for a Ts
%% Create a batch of KnownProfile
% N_SOC = 101; % This will be weird to code
% SOC_min = 0;
% SOC_max = 100;
% 
% SOC_vec = linspace(SOC_min,SOC_max,N_SOC);
% % SOC_vec = SOC_min:5:SOC_max;
% SOC_vec = [54];
% 
% for i = 1:length(SOC_vec)
% %     i
%     % KBSOC    = SOC_vec(i);
%     PRBS_SOC = SOC_vec(i);
%     % EIS_SOC  = SOC_vec(i);
%     % SS_SOC   = SOC_vec(i);
%     CreateProject
% end

%% Loop over Sampling Rates
% N_t_s = 25;
% T_s_min = -1; % T_s = 10^(T_s_min)
% T_s_max =  1; % T_s = 10^(T_s_max)
% 
% T_s_vec = logspace(T_s_min,T_s_max,N_t_s);
% % T_s_vec = T_s_vec(2);
% 
% for i = 1:length(T_s_vec)
%     % i
%     PRBS_Tswitch = 10*T_s_vec(i);
%     CreateProject
% end

%% Soret and Beta loops
    % % ST_vec   = [-1.5];
    % % Beta_vec = [-1]*1e-3;
    % % ITG_vec  = [-2];
    % 
    % % ST_vec   = [-1.5, -1, 0, 1, 1.5];
    % % Beta_vec = [-1, -0.5, 0, 0.5, 1]*1e-3;
    % % ITG_vec  = [-6,-4,-2,2,4,6];          % ITG = T_AN - T_CA
    % 
    % % ST_vec   = [1.5, 2.0, 2.5, 3.0, 3.5, 4];
    % % Beta_vec = [-3, -2.5, -2, -1.5, -1]*1e-3;
    % % ITG_vec  = [-6,-4,-2,2,4,6];          % ITG = T_AN - T_CA
    % 
    % ST_vec   = [4];
    % Beta_vec = [1]*1e-3;
    % ITG_vec  = [-6,-4,-2,2,4,6];          % ITG = T_AN - T_CA
    % 
    % for ss = 1:length(ST_vec)
    %     for bb = 1:length(Beta_vec)
    %         for ii = 1:length(ITG_vec)
    %             preST   = ST_vec(ss);
    %             preBeta = Beta_vec(bb);
    %             preITG  = ITG_vec(ii);
    % 
    %             battery_name = ['ST' num2str(preST) '_Beta' num2str(preBeta ) '_ITG' num2str(preITG)];
    % 
    %             ITGAbsHalf = abs(preITG/2);
    %             if preITG<0
    %                 SIM.preAN_Temp = 20 + 273.15 - ITGAbsHalf;
    %                 SIM.preCA_Temp = 20 + 273.15 + ITGAbsHalf;
    %             else
    %                 SIM.preAN_Temp = 20 + 273.15 + ITGAbsHalf;
    %                 SIM.preCA_Temp = 20 + 273.15 - ITGAbsHalf;
    %             end
    %             SIM.preS_T  = preST;
    %             SIM.preBeta = preBeta;
    %             SIM.preITG  = preITG;
    % 
    %             CreateProject
    %         end
    %     end
    % end

%% Constant Deg Cycling Loops
    % % DegType_vec     = {'ConstantAML'};
    % % DegLocation_vec = {'AN'};
    % % N_cycles_vec    = [4];
    % 
    % DegType_vec     = {'ConstantAML','CombinedCTRG_AML'};
    % DegLocation_vec = {'AN', 'CA', 'Both'};
    % N_cycles_vec    = [4, 20];
    % 
    % % DegType_vec     = {'ConstantCTRG','ConstantAML','CombinedCTRG_AML'};
    % % DegLocation_vec = {'AN', 'CA', 'Both'};
    % % N_cycles_vec    = [4, 20];
    % 
    % for tt = 1:length(DegType_vec)
    %     for kk = 1:length(DegLocation_vec)
    %         for cc = 1:length(N_cycles_vec)
    %             preDegType     = DegType_vec{tt};
    %             preDegLocation = DegLocation_vec{kk};
    %             preN_cycles    = N_cycles_vec(cc);
    % 
    %             battery_name = ['Iso20_' preDegType '_' preDegLocation];
    %             KBCPProfileFilename = ['Chg_Dchg_' num2str(preN_cycles) 'Cycles'];
    % 
    %             if     strcmp(preDegType , 'ConstantCTRG')
    %                 SIM.preCTRGrowth = 1;
    %                 SIM.preAMLoss    = 0;
    %             elseif strcmp(preDegType , 'ConstantAML')
    %                 SIM.preCTRGrowth = 0;
    %                 SIM.preAMLoss    = 1;
    %             elseif strcmp(preDegType , 'CombinedCTRG_AML') 
    %                 SIM.preCTRGrowth = 1;
    %                 SIM.preAMLoss    = 1;
    %             end
    % 
    %             if     strcmp(preDegLocation , 'AN')
    %                 SIM.preCTRG_AN = 1;
    %                 SIM.preCTRG_CA = 0;
    %                 SIM.preAML_AN  = 1;
    %                 SIM.preAML_CA  = 0;
    %             elseif strcmp(preDegLocation , 'CA')
    %                 SIM.preCTRG_AN = 0;
    %                 SIM.preCTRG_CA = 1;
    %                 SIM.preAML_AN  = 0;
    %                 SIM.preAML_CA  = 1;
    %             elseif strcmp(preDegLocation , 'Both') 
    %                 SIM.preCTRG_AN = 1;
    %                 SIM.preCTRG_CA = 1;
    %                 SIM.preAML_AN  = 1;
    %                 SIM.preAML_CA  = 1;
    %             end
    % 
    %             CreateProject
    %         end
    %     end
    % end


%% 0-Crate Tests
    SOC_vec = [10];
    % DeltaT  = -6;
    % SOC_vec = [10, 25, 50, 75, 90];
    DeltaT  = -6:2:6;
    
    for kk = 1:length(SOC_vec)
        for tt = 1:length(DeltaT)
            SOC_start = SOC_vec(kk);

            battery_name = ['dVdT_DelT' num2str(DeltaT(tt)) '_startSOC' num2str(SOC_start)];
            SIM.FLAG_TempBC = 0;

            ITGAbsHalf = abs(DeltaT(tt)/2);
            if DeltaT(tt)<0
                SIM.preAN_Temp = 20 + 273.15 - ITGAbsHalf;
                SIM.preCA_Temp = 20 + 273.15 + ITGAbsHalf;
            else
                SIM.preAN_Temp = 20 + 273.15 + ITGAbsHalf;
                SIM.preCA_Temp = 20 + 273.15 - ITGAbsHalf;
            end
            
            CreateProject
        end
    end
    

