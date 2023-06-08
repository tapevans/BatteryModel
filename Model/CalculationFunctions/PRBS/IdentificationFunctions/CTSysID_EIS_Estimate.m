function [xout,np,nz]=CTSysID_EIS_Estimate(sys,info_all,smoothdat,Iterations)
%CTSYSID_EIS_ESTIMATE
%
% Estimate EIS from multiple experiments using parameterized continuous time model
%
% [xout,np,nz]=CTSysID_EIS_Estimate(sys,info_all,smoothdat)
%
% inputs: 
%   sys - length k cell array of initial system estimates
%   info_all - length k cell array with elements
%        .delayCount - number of samples that PRBS is held constant
%   smoothdata - lengtk k cell array of iddata elements with pre-filtered current/voltage responses

for i=1:length(sys)
    Update_Freq(i) = 1./(info_all{i}.delayCount*sys{i}.Ts);
end

[~,index]=sort(Update_Freq,'ascend');

 [zt,pt,k1]=zpkdata(tf(d2c(sys{index(1)})));
 zt=zt{1};
 [~,i]=sort(real(zt),'descend');
 zt = zt(i);
% z = real(zt(abs(mod(angle(zt),2*pi)-pi)<.2));
 pt=pt{1};
 [~,i]=sort(real(pt),'descend');
 pt=pt(i);

 ind = find(abs(mod(angle(pt),2*pi)-pi)<.2 | abs(pt)<1e-5);
 
 p1 = real(pt(ind));
 z1 = real(zt(ind));
 %k = bode(sys{index(1)},0)*prod(abs(p))/prod(abs(z));
 z=z1;
 p=p1;
 k=k1;
 for i=2:length(sys)
   [zt2,pt2]=zpkdata(tf(d2c(sys{index(i)})));
   zt2=zt2{1};
   [~,i]=sort(real(zt2),'descend');
   zt2=zt2(i);
   pt2=pt2{1};
   [~,i]=sort(real(pt2),'descend');
   pt2=pt2(i);
 
   inew = find(abs(mod(angle(pt2),2*pi)-pi)<.2 & abs(pt2)>.001/sys{index(end)}.Ts & real(pt2)<min(p));
   znew=real(zt2(inew));
   pnew=real(pt2(inew));
   %
   % check that pole magnitude is less than zero magnitude. If not, place zero at 5% greater magnitude
   %
   for i=1:length(pnew)
       if abs(pnew(i))>abs(znew(i)),
           znew(i) = pnew(i)*1.05;
       end;
   end;
           
 z = [z;znew];
 p = [p;pnew];
 k = k*prod(abs(pnew))/prod(abs(znew));
 end;
% x  = [k*prod(pnew)/prod(znew);sort(p);-sort(p)+sort(z)]
 x  = [k;sort(p);sort(z)];

np=length(p);
nz=length(z);

lb = [-inf;-inf*ones(np,1);-inf*ones(nz,1)];
ub =  [inf;0*ones(np,1);-0*ones(nz,1)];
fun =@(x)pzfit_tot(x,np,nz,smoothdat,sys);
options = optimoptions('fmincon');
options.Display = 'iter';
options.Algorithm = 'sqp';
options.OptimalityTolerance = 1e-6;
options.MaxFunctionEvaluations=Iterations;
figure
[xout,fval,exitflag,output,lambda,grad,hessian]=fmincon(fun,x,[],[],[],[],lb,ub,[],options);