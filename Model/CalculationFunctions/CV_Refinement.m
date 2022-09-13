%% Refine Current Profile for Plating Constraint
% Description: 
% 
% Inputs:
%   *
%   *
%   *
%
% Outputs:
%   *
%   *
%%
function [continueRefinement, current_vec , current , sum_of_refinement] = CV_Refinement(cell_voltage, t_soln, SIM, N)
current = SIM.region_current_vec;
continueRefinement = 1;

% Time and associated region
    time_step_regions = zeros(length(t_soln),1);
    region = 1;
    for i = 1:length(t_soln)
        if t_soln(i) < SIM.region_time_vec(region+1)
            time_step_regions(i) = region;
        else
            region = region + 1;
            time_step_regions(i) = region;
        end
    end
    
% Refinement
    no_refine = zeros(SIM.N_regions,1);
    for r = 1:SIM.N_regions
        if r > 90
            cell_voltage; % Troubleshooting
        end
        idx = find(time_step_regions == r);
        min_del_phi = min(cell_voltage(idx,N.N_CV_AN));
        if min_del_phi < 0 % Delta Phi constraint is violated
            current(r) = SIM.region_current_vec(r) - SIM.i_user_OG * SIM.decrease_percent / 100;
            if abs(current(r)) < abs(SIM.C_rate_min * SIM.i_user_OG) % Saturation limit
                current(r) = SIM.C_rate_min * SIM.i_user_OG;
            end
        elseif r > max(time_step_regions)
            no_refine(r) = 1;
            current(r) = current(r-1);
        else
            if min_del_phi < SIM.tol_Delta_phi % No adjustment needed, within tolerance
                no_refine(r) = 1;
            else % Increase C-rate to get closer to plating curve
                current(r) = SIM.region_current_vec(r) + SIM.i_user_OG * SIM.increase_percent / 100;
            end
        end
    end
       
    % Remake the current vector with the refinement
    current_vec = 0;

    for i = 1:SIM.N_regions-1
        current_vec(end+1) = current(i);
        current_vec(end+1) = current(i);
        current_vec(end+1) = (current(i) + current(i+1))/2;
    end

    i = SIM.N_regions;
    current_vec(end+1) = current(i);
    current_vec(end+1) = current(i);
    
% More Refinement Needed
    sum_of_refinement = sum(no_refine);
    if sum_of_refinement == SIM.N_regions
        continueRefinement = 0;
    end
    
end