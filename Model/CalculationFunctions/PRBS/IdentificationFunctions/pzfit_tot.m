function error = pzfit_tot(x,np,nz,data,sys)
k = x(1);
p = x(1+(1:np));
z = x(1+np+(1:nz));
if length(z)<length(p),
    z = [z;zeros(length(p)-length(z),1)];
end;
G = minreal(zpk(z,p,k));
clf
error=0;

for i=1:length(data);
 Ts = data{i}.Ts;
 u = data{i}.u;
 y = data{i}.y;
 t = Ts*(0:length(u)-1);
 if (1),
  % no noise case
  yhat = lsim(G,u,t,[],'zoh');
%  error = error + norm(y - yhat)^2/(sys{i}.NoiseVariance*length(u));
%  error = error + norm(y - yhat)^2/(10e-10*length(u));
  error = error + norm(y - yhat)^2/(10e-10*length(u)*max(abs(y)));
 else
   % noise case
 [e,yhat]=KFerrorfit(G,sys{i}.Ts,y,u,1,1);
 error = error + e/(sys{i}.NoiseVariance*length(u));
 end;
 subplot(length(data)+1,1,i);
 plot(t,y,t,yhat)
end;
subplot(length(data)+1,1,length(data)+1)
semilogx(-p,0*p,'x')
hold on
semilogx(-z,0*p,'o')
hold off
drawnow


A = [zeros(np,1),-eye(np),eye(np)];

b = zeros(np,1);
%disp(['constraint cost:',num2str(0*sum(max(A*x-b,0)))])
error = error + 10*sum(max(A*x-b,0));
