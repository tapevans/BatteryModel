%% getImpedanceFromSSSystem

function [Z_results] = getImpedanceFromSSSystemAllOutputs_pinv(sys , M , omega , SIM , P)
%% Initialize
    A = sys.A;
    B = sys.B;
    C = sys.C;
    D = sys.D;

    [N_out,~] = size(C);

    Z_all = cell(N_out,1);
    Z = zeros(length(omega),1);


%% Initialize Cell Array
    for j = 1:N_out
        Z_all{j} = Z;
    end


 %% Calculate Impedance   
    for i = 1:length(omega)
        s      = omega(i)*(1i);
        Z_temp = C * (pinv(s*M - A)*B) + D;
        for j = 1:N_out
            Z_all{j}(i) = Z_temp(j);
        end
        % Z(i)   = Z_temp(1); % Hard-coded for cell voltage
    end
    multiple = -SIM.A_c^-1;
    for j = 1:N_out
        Z_new_all{j,1} = multiple*Z_all{j};
    end


%% Combine results to a single variable
for j = 1:N_out
    Z_new = Z_new_all{j,1};
% Calculations
    Z_mag = sqrt(real(Z_new).^2 + imag(Z_new).^2);
    Z_dB = 20*log10(Z_mag);
    Z_angle_deg = asind(imag(Z_new)./Z_mag);
    
% Save
    Z_results{j,1}(:,P.SS.omega)    = omega';
    Z_results{j,1}(:,P.SS.Z_mag)    = Z_mag;
    Z_results{j,1}(:,P.SS.Z_Re)     = real(Z_new);
    Z_results{j,1}(:,P.SS.Z_Im)     = imag(Z_new);
    Z_results{j,1}(:,P.SS.Z_dB)     = Z_dB;
    Z_results{j,1}(:,P.SS.Z_ps_deg) = Z_angle_deg;

    % if Z_results(1,P.SS.Z_Re) < 0 
    %     Z_results(:,P.SS.Z_Re) = -1 * Z_results(:,P.SS.Z_Re);
    % end
end
end