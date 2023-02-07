function [sys_CT3 , sys_DT3 , sys_CT5 , sys_DT5] = getSS_System(SIM,N,P,FLAG)
% 3 State System 
    [A,B,C,D] = getAll_SS3(SIM,N,P,FLAG);
    
    sys = dss(A,B,C,D,SIM.Mass3);
    sys_CT3 = ss(sys,'explicit');
    sys_DT3 = c2d(sys_CT3,SIM.Ts);

% 5 State System
    [A,B,C,D] = getAll_SS5(SIM,N,P,FLAG);
    
    sys = dss(A,B,C,D,SIM.Mass5);
    sys_CT5 = ss(sys,'explicit');
    sys_DT5 = c2d(sys_CT5,SIM.Ts);

end