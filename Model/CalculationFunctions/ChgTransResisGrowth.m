%% Charge Transfer Resistance Growth
% This function calculates the growth rate of the charge transfer
% resistance. 
function [CTR_GrowthRate] = ChgTransResisGrowth(del_phi , SIM__delPhi_threshold , SIM__CTRG_cons_vec, FLAG__CTRG_Dep)
    if FLAG__CTRG_Dep
        del_phi_Condition = del_phi < SIM__delPhi_threshold;
        CTR_GrowthRate    = SIM__CTRG_cons_vec .* abs(del_phi-SIM__delPhi_threshold) .* del_phi_Condition;
        
        % R_SEI_GrowthRate = [AN__ChgTranResGrowthRate*ones(1,N.N_CV_AN) , zeros(1 , N.N_CV_SEP ) , CA__ChgTranResGrowthRate*ones(1,N.N_CV_CA)];
        % R_SEI_GrowthRate = R_SEI_GrowthRate .* abs(del_phi-delPhi_threshold); %% Could just throw a negative sign in front of this instead of using abs()
    else
        CTR_GrowthRate    = SIM__CTRG_cons_vec;
    end
end

    % del_phi = SV( P.del_phi,: );
    % delPhi_threshold = 0.1; %%%%%!!!!!!!!!! Hardcoded
    % del_phi_neg = del_phi < delPhi_threshold;
    % R_SEI_GrowthRate = [AN.R_SEI_growthRate*ones(1,N.N_CV_AN) , zeros(1 , N.N_CV_SEP ) , CA.R_SEI_growthRate*ones(1,N.N_CV_CA)];
    % R_SEI_GrowthRate = R_SEI_GrowthRate .* abs(del_phi-delPhi_threshold); %% Could just throw a negative sign in front of this instead of using abs()
    % R_SEI_GrowthRate = R_SEI_GrowthRate .* del_phi_neg * FLAG.SEI_growth;