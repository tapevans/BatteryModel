%% Batt Residual Function
% This function is used to return the time derivative of the governing equations.

function dSVdt = batt_GovEqn_test(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user)
%% Organize (reshape) the SV
SV = SV1Dto2D(SV , N , P , FLAG);
SV = addPhiEl2SV(SV,P,N);

%% Obtain Property Values
if FLAG.VARIABLE_PROPS_FROM_HANDLES
    props = getProps( SV , AN , SEP, CA , EL , P , N , CONS , FLAG , PROPS);
else
    props = PROPS;
end

%% Calculate i_user
if SIM.SimMode == 3 % State Space EIS
    % Uses the i_user value from the function handle
elseif SIM.SimMode == 4 % KBCP Controller
    % Uses the i_user value from the function handle
elseif SIM.SimMode == 5 % MOO Controller
    % Uses the i_user value from the function handle
else
    i_user = i_user_calc(t,SIM);
end

%% Calculate All Fluxes
[i_ed , i_el ] = currentCalc( SV , AN , SEP , CA , EL , P , N , CONS , FLAG , i_user , props);

i_Far = iFarCalc( SV , AN , CA , P , N , CONS , FLAG , props );
s_dot = i_Far / CONS.F ;

J_Liion = JLiionCalc( SV , AN , SEP, CA , EL , P , N , CONS , FLAG , i_el , props);
J_Li    = JLiCalc( SV , AN , CA , P , N , s_dot , props);

%%
q_cond = zeros(1 , N.N_CV_tot+1);
q_conv = zeros(1 , N.N_CV_tot  );
q_gen  = zeros(1 , N.N_CV_tot  );
if FLAG.COE
    % CALCULATE CONDUCTION WITH BC
    q_cond = qCondCalc( SV , AN , SEP, CA , EL , P , N , SIM , CONS , FLAG , props);
    
    if FLAG.CV_CONV % CALCULATE SURFACE CONVECTION
        %q_conv 
    end
    
    if FLAG.HEAT_GEN_TOTAL % TOTAL HEAT GEN CALC
        q_gen = calcHeatGenTot( SV , AN , SEP, CA , EL , P , N , CONS , SIM , FLAG , props , i_el, i_ed, i_Far);
    end
end

%% Calculate Equilibrium Voltage
E_eq_an  = AN.EqPotentialHandle( SV(P.C_Li_surf_AN , N.CV_Region_AN ) / AN.C_Li_max );
E_eq_ca  = CA.EqPotentialHandle( SV(P.C_Li_surf_CA , N.CV_Region_CA ) / CA.C_Li_max );

E_eq_vec = [ E_eq_an , zeros(1,N.N_CV_SEP) , E_eq_ca];


%% Solving for dSVdt
%% ---- Anode ----
dSVdt_AN = zeros(N.N_SV_AN,N.N_CV_AN);
for i = 1:N.N_CV_AN
    % Temp
%     dSVdt_an(P.T, i) = (  -(q_cond(i+1) - q_cond(i))/ AN.del_x...
%                           - q_conv(i)...
%                           + q_gen(i)  )...
%                         / ( props(P.rho,i) * props(P.c_p,i) );
    dSVdt_AN(P.T, i) = -(SV(P.T,i) - SIM.T_inf);
    
    % phi_el                    
    dSVdt_AN(P.del_phi , i) =  (AN.A_c / AN.A_surf_CV)*(i_el(i+1) - i_el(i)  ) ...
                                - (SV(P.V_1,i) - SV(P.phi_el,i))/AN.R_SEI;
                        
    % phi_ed                                        
    dSVdt_AN(P.phi_ed , i) =   i_ed(i  ) + i_el(i  )...
                             - (i_ed(i+1) + i_el(i+1));
    
    % V_1                    
    dSVdt_AN(P.V_1    , i) = - (SV(P.V_1,i) - SV(P.phi_el,i))/AN.R_SEI ...
                             +  i_Far(i);
                        
    % V_2                    
    dSVdt_AN(P.V_2    , i) =  - i_Far(i) ...
                              + SV(P.i_PS,i);
    
	% i_PS                    
    dSVdt_AN(P.i_PS   , i) =    SV(P.phi_ed,i) - SV(P.V_2,i) -  E_eq_vec(i);
                        
	% C_Li^+
    dSVdt_AN(P.C_Liion, i) = -(J_Liion(i+1) - J_Liion(i))/ AN.del_x...
                              + s_dot(i) * AN.A_s;
                
    % C_Li
    for j = 1:N.N_R_AN
            dSVdt_AN(P.C_Li+j-1, i) = -3*(AN.r_half_vec(j+1)^2 * J_Li(j+1,i) - AN.r_half_vec(j)^2 * J_Li(j,i)) ...
                                       / (AN.r_half_vec(j+1)^3-AN.r_half_vec(j)^3);
    end
end
% Fix Boundary Conditions
i = 1;
dSVdt_AN(P.phi_ed, i) = -SV(P.phi_ed,i);

%% ---- Separator ----
dSVdt_SEP = zeros(N.N_SV_SEP,N.N_CV_SEP);
for i = 1:N.N_CV_SEP 
    index_offset = N.N_CV_AN + i;  
    % Temp
    dSVdt_SEP(P.SEP.T, i) = -(SV(P.T,index_offset) - SIM.T_inf);
%     dSVdt_SEP(P.SEP.T, i) = (  -(q_cond(index_offset+1) - q_cond(index_offset))/ SEP.del_x...
%                                - q_conv(index_offset)...
%                                + q_gen(index_offset)  )...
%                             / ( props(P.rho,index_offset) * props(P.c_p,index_offset) );
    
    % phi_el
    dSVdt_SEP(P.SEP.phi_el, i) = -(i_el(index_offset+1)-i_el(index_offset));
    
    % C_Li^+
    dSVdt_SEP(P.SEP.C_Liion, i)= -(J_Liion(index_offset+1) - J_Liion(index_offset))/SEP.del_x;
end

%% ---- Cathode ----
dSVdt_CA = zeros(N.N_SV_CA,N.N_CV_CA);
for i = 1:N.N_CV_CA
    index_offset = N.N_CV_AN + N.N_CV_SEP + i;  
    % Temp
    dSVdt_CA(P.T, i) =  -(SV(P.T,index_offset) - SIM.T_inf);
% 	dSVdt_CA(P.T, i) =  (  -(q_cond(index_offset+1) - q_cond(index_offset))/ CA.del_x...
%                            - q_conv(index_offset)...
%                            + q_gen(index_offset)  )...
%                            / ( props(P.rho,index_offset) * props(P.c_p,index_offset) );
    
    % phi_el
    dSVdt_CA(P.del_phi , i) =  (CA.A_c / CA.A_surf_CV)*(i_el(index_offset+1) - i_el(index_offset)  ) ...
                                - (SV(P.V_1,index_offset) - SV(P.phi_el,index_offset))/CA.R_SEI;
    
	% phi_ed                    
    dSVdt_CA(P.phi_ed , i)  =   i_ed(index_offset  ) + i_el(index_offset  )...
                              - (i_ed(index_offset+1) + i_el(index_offset+1));
    
	% V_1                    
    dSVdt_CA(P.V_1    , i)  = - (SV(P.V_1,index_offset) - SV(P.phi_el,index_offset))/CA.R_SEI ...
                              +  i_Far(index_offset);
    
	% V_2                    
    dSVdt_CA(P.V_2    , i)  = -  i_Far(index_offset)...
                              +  SV(P.i_PS,index_offset);
    
	% i_PS                    
    dSVdt_CA(P.i_PS   , i)  =    SV(P.phi_ed,index_offset) - SV(P.V_2,index_offset) ...
                              -  E_eq_vec(index_offset);
    
    % C_Li^+
    dSVdt_CA(P.C_Liion, i) = -(J_Liion(index_offset+1) - J_Liion(index_offset))/ CA.del_x...
                              + s_dot(index_offset) * CA.A_s;
    % C_Li
    for j = 1:N.N_R_CA
            dSVdt_CA(P.C_Li+j-1, i) = -3*(CA.r_half_vec(j+1)^2 * J_Li(j+1,index_offset) - CA.r_half_vec(j)^2 * J_Li(j,index_offset)) ...
                                       / (CA.r_half_vec(j+1)^3-CA.r_half_vec(j)^3);
    end
end

%% Reshape
% Reshape matrix to a column vector
dSVdt_AN  = reshape(dSVdt_AN ,[],1);
dSVdt_SEP = reshape(dSVdt_SEP,[],1);
dSVdt_CA  = reshape(dSVdt_CA ,[],1);

%% Combine all dSVdt from each region
dSVdt = [dSVdt_AN ; dSVdt_SEP ; dSVdt_CA];

%% Used for troubleshooting
% if t>SIM.initial_offset
%    t;
% end
% if t>SIM.t_ramp
%    t;
% end
end