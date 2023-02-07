function [InputSignal] = getInputSignal(SIM,N,P,FLAG)
[V_in] = V_inCalc(SIM.t_vec,SIM,FLAG);
InputSignal = [SIM.t_vec' , V_in'];


end