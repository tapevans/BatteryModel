%% Update Properties
% This function will be used to determine the desired properties as a
% function of temperature and concentration
function props = getProps(  Ce, T, X_AN,  X_CA, FLAG__VARIABLE_kappa                                                                                              , ...
                            FLAG__VARIABLE_D_Liion, FLAG__VARIABLE_activity, FLAG__VARIABLE_tf_num, FLAG__VARIABLE_D_o_AN, FLAG__VARIABLE_D_o_CA, FLAG__Bruggeman , ...
                            P__kappa, P__D_o_Li_ion, P__activity, P__tf_num, P__D_o                                                                               , ...
                            N__N_R_AN, N__N_R_CA, N__CV_Region_AN, N__CV_Region_CA                                                                                , ...
                            CONS__BRUG                                                                                                                            , ...
                            AN__D_oHandle, CA__D_oHandle                                                                                                          , ...
                            PROPS, EL__kappaHandle, EL__D_o_Li_ionHandle, EL__ActivityHandle, EL__tf_numHandle                                                         )

%% Initialize       
    props = PROPS;
    
%% Evaluate Properties 
    % kappa
        if FLAG__VARIABLE_kappa
            props( P__kappa , : )      = EL__kappaHandle(Ce , T);
        else
            props( P__kappa , : )      = PROPS( P__kappa , :);
        end
    
    % D_o_Liion
        if FLAG__VARIABLE_D_Liion
            props( P__D_o_Li_ion , : ) = EL__D_o_Li_ionHandle(Ce , T);
        else
            props( P__D_o_Li_ion , : ) = PROPS( P__D_o_Li_ion , :);
        end
    
    % Activity
        if FLAG__VARIABLE_activity
            props( P__activity , : )   = EL__ActivityHandle(Ce , T);
        else
            props( P__activity , : )   = PROPS( P__activity , :);
        end
    
    % Transference Number
        if FLAG__VARIABLE_tf_num
            props( P__tf_num , : )     = EL__tf_numHandle(Ce , T);
        else
            props( P__tf_num , : )     = PROPS( P__tf_num , :);
        end
    
    % Active Material Diffusion -- Anode
        if FLAG__VARIABLE_D_o_AN
            props( P__D_o:(P__D_o+N__N_R_AN-1) , N__CV_Region_AN ) = AN__D_oHandle(X_AN , T(:,N__CV_Region_AN));
        else
            props( P__D_o:end , N__CV_Region_AN ) = PROPS( P__D_o:end , N__CV_Region_AN );
        end
    
    % Active Material Diffusion -- Cathode
        if FLAG__VARIABLE_D_o_CA
            props( P__D_o:(P__D_o+N__N_R_CA-1) , N__CV_Region_CA ) = CA__D_oHandle(X_CA , T(:,N__CV_Region_CA));
        else
            props( P__D_o:end , N__CV_Region_CA ) = PROPS( P__D_o:end , N__CV_Region_CA);
        end
    
    % % Parameters that don't have handles yet %%%%%%%% Can I improve this?
    %     props( P__sigma      , N__CV_Region_AN ) = AN__sigma * ones( 1 , N__N_CV_AN );
    %     props( P__sigma      , N__CV_Region_CA ) = CA__sigma * ones( 1 , N__N_CV_CA );
    %     props( P__R_SEI      , N__CV_Region_AN ) = AN__R_SEI * ones( 1 , N__N_CV_AN );
    %     props( P__R_SEI      , N__CV_Region_CA ) = CA__R_SEI * ones( 1 , N__N_CV_CA );
    %     props( P__lambda_eff , : )              = [AN__lambda_eff * ones(1,N__N_CV_AN) , SEP__lambda_eff * ones(1,N__N_CV_SEP), CA__lambda_eff * ones(1,N__N_CV_CA)];
    %     props( P__rho_eff    , : )              = [AN__rho_eff    * ones(1,N__N_CV_AN) , SEP__rho_eff    * ones(1,N__N_CV_SEP), CA__rho_eff    * ones(1,N__N_CV_CA)];
    %     props( P__c_p_eff    , : )              = [AN__c_p_eff    * ones(1,N__N_CV_AN) , SEP__c_p_eff    * ones(1,N__N_CV_SEP), CA__c_p_eff    * ones(1,N__N_CV_CA)];
        
    % Account for tortuosity (Bruggeman)
        if FLAG__Bruggeman
            props = props .* CONS__BRUG;
        end
end