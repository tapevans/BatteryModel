%% Active Material Loss
% This function calculates the loss rate of the number of particles
function [AM_LossRate] = NParticlesLoss(C_Li_diff , SIM__C_Li_diff_threshold , SIM__AML_cons_vec , FLAG__AML_Dep)
    if FLAG__AML_Dep
        C_Li_Condition = C_Li_diff < SIM__C_Li_diff_threshold;
        AM_LossRate    = SIM__AML_cons_vec .* abs(C_Li_diff-SIM__C_Li_diff_threshold) .* C_Li_Condition;
    else
        AM_LossRate    = SIM__AML_cons_vec;
    end
end