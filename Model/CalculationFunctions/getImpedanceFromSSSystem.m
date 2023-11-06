%% getImpedanceFromSSSystem

function [Z_results] = getImpedanceFromSSSystem(sys , M , omega , SIM , P)
%% Initialize
    A = sys.A;
    B = sys.B;
    C = sys.C;
    D = sys.D;

    Z = zeros(length(omega),1);

 %% Calculate Impedance   
    for i = 1:length(omega)
        s      = omega(i)*(1i);
        Z_temp = C * ((s*M - A)\B) + D;
        Z(i)   = Z_temp(1); % Hard-coded for cell voltage
    end
    multiple = -SIM.A_c^-1;
    Z_new = multiple*Z;


%% Combine results to a single variable
% Calculations
    Z_mag = sqrt(real(Z_new).^2 + imag(Z_new).^2);
    Z_dB = 20*log10(Z_mag);
    Z_angle_deg = asind(imag(Z_new)./Z_mag);
    
% Save
    Z_results(:,P.SS.omega)    = omega';
    Z_results(:,P.SS.Z_mag)    = Z_mag;
    Z_results(:,P.SS.Z_Re)     = real(Z_new);
    Z_results(:,P.SS.Z_Im)     = imag(Z_new);
    Z_results(:,P.SS.Z_dB)     = Z_dB;
    Z_results(:,P.SS.Z_ps_deg) = Z_angle_deg;

    % if Z_results(1,P.SS.Z_Re) < 0 
    %     Z_results(:,P.SS.Z_Re) = -1 * Z_results(:,P.SS.Z_Re);
    % end
end