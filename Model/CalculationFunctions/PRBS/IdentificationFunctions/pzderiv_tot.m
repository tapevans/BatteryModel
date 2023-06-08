function [Gnom,P] = pzderiv_tot(x,np,nz,data,sys)
k = x(1);
p = x(1+(1:np));
z = x(1+np+(1:nz));
Gnom = zpk(z,p,k);
dyhat_dtheta = [];
yhat=[];
for i=1:length(data)
  Ts = data{i}.Ts;
   u = data{i}.u;
   t = Ts*(0:length(u)-1);
  yhat = [yhat;lsim(Gnom,u,t)/sqrt(sys{i}.NoiseVariance)];
end;
for ell=1:length(x)
    del = zeros(size(x));
    del(ell) = max(abs(x(ell))*1e-6,1e-9);
    dx = x-del;
    dk = dx(1);
    dp = dx(1+(1:np));
    dz = dx(1+np+(1:nz));
    Gdel=zpk(dz,dp,dk);
    dyhat=[];
    for i=1:length(data)
       Ts = data{i}.Ts;
       u = data{i}.u;
       t = Ts*(0:length(u)-1);
       dyhat = [dyhat;lsim(Gdel,u,t)/sqrt(sys{i}.NoiseVariance)];
    end;
    dyhat_dtheta(:,ell) = (yhat-dyhat)/del(ell);
end;
%P = lambda*inv(dyhat_dtheta'*dyhat_dtheta);
P = inv(dyhat_dtheta'*dyhat_dtheta);
    
