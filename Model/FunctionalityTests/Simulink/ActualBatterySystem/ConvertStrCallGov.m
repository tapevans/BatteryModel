%% Convert String to Handle
% And then call governing equations function
function [dSVdt] = ConvertStrCallGov(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user)
AN.EqPotentialHandle = str2func(AN.EqPotentialHandle);
AN.i_oHandle         = str2func(AN.i_oHandle);
AN.sigmaHandle       = str2func(AN.sigmaHandle);
AN.D_oHandle         = str2func(AN.D_oHandle);

CA.EqPotentialHandle = str2func(CA.EqPotentialHandle);
CA.i_oHandle         = str2func(CA.i_oHandle);
CA.sigmaHandle       = str2func(CA.sigmaHandle);
CA.D_oHandle         = str2func(CA.D_oHandle);

EL.tf_numHandle      = str2func(EL.tf_numHandle);
EL.ActivityHandle    = str2func(EL.ActivityHandle);
EL.D_o_Li_ionHandle  = str2func(EL.D_o_Li_ionHandle);
EL.kappaHandle       = str2func(EL.kappaHandle);

dSVdt = batt_GovEqn_test(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user);

end