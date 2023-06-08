function [sysc,sys,totdat,iddat]=IdentificationProcess(V,C,t,Update_Freq,dofilter)
%IDENTIFICATIONPROCESS
%
%  [sysc,sys,totdat,iddat,smoothdat]=IdentificationProcess(V,C,t,Update_Freq) 
%
% inputs:
%    V - N x 1 vector of voltage measurements
%    C - N x 1 vector of current measurements
%    t - N x 1 vector of time samples. It is assumed that the sample time (Ts) is fixed and equal to t(2)-t(1)
%    Update_Freq - 1 / sample time of PRBS sequence  (i.e. smallest time for constant PRBS signal)
%
% outputs: 
%    sysc - estimate as continuous time system (inverse zero order hold)
%    sys - estimate as discrete time system 
%    totdat - iddata structure of data, filtered by pre-filter. Pre-filter band is pi*Update_Freq/1000 to pi/Ts
%    idata - iddata structure of data used for identification of cross validation (first 2/3rds)
%
Ts=mean(diff(t));
 
filtband = [pi*Update_Freq/1000 pi*Update_Freq];
N = length(V);

Nzero = min(find(C~=0))-1;
if Nzero ==0, Nzero=1; end;
Vnom = mean(V(1:Nzero));
y=V-Vnom;
u=C;

totdat = iddata(y,u,Ts);
if dofilter,
    totdat = idfilt(totdat,filtband);
end;
%
% Last 1/3 of data for validation
%
valdat = iddata(y(ceil(2*N/3):N),u(ceil(2*N/3):N),Ts);
if dofilter,
    valdat = idfilt(valdat,filtband);
end;
%
% First 2/3 of data for identification
%
iddat = iddata(y(1:floor(2*N/3)),u(1:floor(2*N/3)),Ts);
arx_opt = arxOptions('InitialCondition','Estimate','EnforceStability',true);
if dofilter,
    iddat = idfilt(iddat,filtband);
end;
%
% Identify models with different orders and use validation data to
% choose correct order
%
m=0;
orderlist=2:8;
for order=orderlist,
    m=m+1;
    arx_sys = arx(iddat,[order order+1 0],arx_opt);
    init_sys=idpoly(1,arx_sys.b,1,1,arx_sys.a,1,iddat.ts);
    oe_opt=oeOptions('InitialCondition','Estimate','EstCovar',true,'EnforceStability',true);
    sys_order{m} = oe(iddat,init_sys,oe_opt);
    if abs(eig(sys_order{m}))>1,
        oe_opt=oeOptions('InitialCondition','Estimate','EstCovar',true,'Focus','Stability');
        sys_order{m} = oe(iddat,init_sys,oe_opt);
    end;
    compare_opt = compareOptions('InitialCondition','e');
    %compare(valdat{k},sys_order{m},compare_opt);
    [yh, fit, X0]= compare(valdat,sys_order{m},compare_opt);
    valfit(m)=fit;
end;
[vmin,imin]=max(valfit);
%
% Set initial model using best model order
%
arx_sys = arx(iddat,[orderlist(imin) orderlist(imin)+1 0],arx_opt);
%if ~measured_current
%    %
%    % Re-compute smoothed voltate estimate and re-compute current
%    %
%    yh = compare(dat,arx_sys,compare_opt);
%    u_smooth=-(Vnom+yh.y).*states(1:floor(2*N/3))/Resistance;
%    smoothdat = iddata(y(1:floor(2*N/3)),u_smooth,mean(diff(t)));
%    if dofilter,
%      smoothdat = idfilt(smoothdat,filtband);
%    end;
%else
%    smoothdat = iddat;
%end;

smoothdat=totdat;


%
% Identify model using smoothed input
%
B = arx_sys.b;
A = arx_sys.a;
init_sys=idpoly(1,B,1,1,A,1,iddat.ts);
oe_opt=oeOptions('InitialCondition','Estimate','EstCovar',true);
%    oe_opt=oeOptions('InitialCondition','zero','EstCovar',true);
sys=oe(smoothdat,init_sys,oe_opt);
if abs(eig(sys))>1,
    oe_opt=oeOptions('InitialCondition','Estimate','EstCovar',true,'Focus','Stability');
    sys{m} = oe(smoothdat,init_sys,oe_opt);
end;
if (0),
    figure
    compare(sys,dat)
end;
%
% convert estimated model to continous time
%
sysc = d2c(sys,'zoh');