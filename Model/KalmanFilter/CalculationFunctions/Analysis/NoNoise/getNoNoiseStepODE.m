function [t_soln, x_soln, z_soln] = getNoNoiseStepODE(SIM,FLAG,N,P)
% Get Filename
    filename = getODEStepFilename(FLAG);
    data = load(filename);
    
% % Plot Data
% z = OutputMatrix*data.SOLN.y;
% figure
% plot(data.SOLN.x,z(1,:),'-o')

% Get just the CC part and then add on the rest later
    t_final_Ts = SIM.Ts*FLAG.N_samples;
    t_CC_vec = 0:SIM.Ts:t_final_Ts;
    
    x = deval(data.SOLN,t_CC_vec);
    
    z = SIM.OutputMatrix * x;

% Add Rest Data
    t_sim_vec = [0 SIM.Ts (t_CC_vec+2*SIM.Ts)];
    x = [x(:,1), x(:,1) , x];
    z = [z(:,1), z(:,1) , z];

% Convert Return Variables
    t_soln = t_sim_vec;
    x_soln = x;
    z_soln = z;

end