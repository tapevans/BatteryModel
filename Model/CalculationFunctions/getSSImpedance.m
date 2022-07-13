%% getSSImpedance
% This function is used to calculate the state space matricies and its
% associated impedance.

%%%%%%%%%%%%%Add more here about what A,B,C,D are and the SS equations

function [A,B,C,D,Z_results] = getSSImpedance(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS)
%% Initializing Inputs
SV      = SIM.SV_IC;
i_user  = SIM.i_user;
M       = SIM.M;
freq    = SIM.freq;
TOL.Rel = 1e-3; % Relative tolerance for perturbation
TOL.Abs = 1e-6; % Absolute tolerance for perturbation

%% Solve for State Space (SS) Matricies
[A,B,C,D] = getSSMatricies(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,SV,i_user,TOL);

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

%% State-Space Matricies
function [A,B,C,D] = getSSMatricies(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,SV,i_user,TOL)
% Initialize
A = zeros(N.N_SV_tot , N.N_SV_tot);
B = zeros(N.N_SV_tot , N.N_In    );
C = zeros(N.N_Out    , N.N_SV_tot);
D = zeros(N.N_Out    , N.N_In    );
t = 0;

inputs_vec = i_user;

% Initial dSVdt
dSVdt_init  = batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec);
output_init = get_output( t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec);

% Perturb SV
for i = 1:N.N_SV_tot
    p    = zeros(N.N_SV_tot,1);
    p(i) = TOL.Rel * SV(i) + TOL.Abs;
    SV_p = SV + p;
    dSVdt = batt_GovEqn(t,SV_p,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec);
        
    A(:,i) = ( dSVdt - dSVdt_init ) / p(i);
    
    output = get_output(t,SV_p,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec);
    C(:,i) = ( output - output_init ) / p(i);
end

% Perturb Inputs
for i = 1:N.N_In
    p    = zeros(N.N_In,1);
    p(i) = TOL.Rel * inputs_vec(i) + TOL.Abs; 
    inputs_vec_p = inputs_vec + p;
    dSVdt = batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec_p);
        
    B(:,i) = ( dSVdt - dSVdt_init ) / p(i);
    
    output = get_output(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec_p);
    D(:,i) = ( output - output_init ) / p(i);
end

end

%% Output
function [output] = get_output(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec)
    output = SIM.OutputMatrix*SV;
end