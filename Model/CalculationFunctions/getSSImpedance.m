%% getSSImpedance
% This function is used to calculate the state space matricies and its
% associated impedance.
function [A,B,C,D,Z_results] = getSSImpedance(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS)
%% Initializing Inputs
SV      = SIM.SV_IC;
i_user  = SIM.i_user;
M       = SIM.M;
freq    = SIM.freq;

%% Solve for State Space (SS) Matricies
[A,B,C,D] = getSSMatricies(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,SV,i_user);

%% Calculate Impedance
Z = zeros(length(freq),1);
for i = 1:length(freq)
    s      = freq(i)*(1i);
    Z_temp = C * ((s*M - A)\B) + D;
    Z(i)   = Z_temp(1); %%%%%%%% Hard-coded for cell voltage
end
multiple = -SIM.A_c^-1;
Z_new = multiple*Z;

%% Combine results to a single variable
    % Calculations
    Z_mag = sqrt(real(Z_new).^2 + imag(Z_new).^2);
    Z_dB = 20*log10(Z_mag);
    Z_angle_deg = asind(imag(Z_new)./Z_mag);
    
    % Save
    Z_results(:,P.SS.omega)    = freq';
    Z_results(:,P.SS.Z_mag)    = Z_mag;
    Z_results(:,P.SS.Z_Re)     = real(Z_new);
    Z_results(:,P.SS.Z_Im)     = imag(Z_new);
    Z_results(:,P.SS.Z_dB)     = Z_dB;
    Z_results(:,P.SS.Z_ps_deg) = Z_angle_deg;

end