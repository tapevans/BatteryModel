%% getSSMatricies
% The purpose of this function is to obtain a state space (SS)
% representation of a lithium-ion battery at a given state. 
%
% Inputs:
% - AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS are the structs containing
%   information about the battery model. 
% - SV is the state vector to be used for the identification
% - i_user is the current flux [A m^-2] used for the identification
%
% Outputs:
% - A,B,C,D are the identified SS matricies

% - Could add more about the dimensions of the matricies
%
%
%
function [A,B,C,D] = getSSMatricies(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,SV,i_user)
%% Parameters
% TOL.Rel = 1e-3; % Relative tolerance for perturbation
% TOL.Abs = 1e-6; % Absolute tolerance for perturbation
TOL.Rel = 1e-5; % Relative tolerance for perturbation
TOL.Abs = 1.5e-8; % Absolute tolerance for perturbation

FLAG.CentralDiff = 1;

%% Initialize
A = zeros(N.N_SV_tot , N.N_SV_tot);
B = zeros(N.N_SV_tot , N.N_In    );
C = zeros(N.N_Out    , N.N_SV_tot);
D = zeros(N.N_Out    , N.N_In    );
t = 0;

inputs_vec = i_user;

%% Initial dSVdt
dSVdt_init  = batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec);
output_init = get_output( t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec);

%% Perturb SV
for i = 1:N.N_SV_tot
    p    = zeros(N.N_SV_tot,1);
    p(i) = TOL.Rel * SV(i) + TOL.Abs;
    SV_p = SV + p;
    dSVdt_p  = batt_GovEqn(t,SV_p,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec);
    output_p = get_output( t,SV_p,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec);
    
    if FLAG.CentralDiff
        SV_m = SV - p;
        dSVdt_m  = batt_GovEqn(t,SV_m,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec);
        output_m = get_output( t,SV_m,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec);
        
        A(:,i) = ( dSVdt_p  - dSVdt_m  ) / (2*p(i));
        C(:,i) = ( output_p - output_m ) / (2*p(i));
    else
        A(:,i) = ( dSVdt_p  - dSVdt_init  ) / p(i);
        C(:,i) = ( output_p - output_init ) / p(i);
    end
end

%% Perturb Inputs
for i = 1:N.N_In
    p    = zeros(N.N_In,1);
    p(i) = TOL.Rel * inputs_vec(i) + TOL.Abs; 
    inputs_vec_p = inputs_vec + p;
    dSVdt_p  = batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec_p);
    output_p = get_output( t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec_p);

    if FLAG.CentralDiff
        inputs_vec_m = inputs_vec - p;
        dSVdt_m  = batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec_m);
        output_m = get_output( t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec_m);

        B(:,i) = ( dSVdt_p  - dSVdt_m  ) / (2*p(i));
        D(:,i) = ( output_p - output_m ) / (2*p(i));
    else
        B(:,i) = ( dSVdt_p  - dSVdt_init  ) / p(i);
        D(:,i) = ( output_p - output_init ) / p(i);
    end
end

end


%% Output
function [output] = get_output(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,inputs_vec)
    output = SIM.OutputMatrix*SV;
end