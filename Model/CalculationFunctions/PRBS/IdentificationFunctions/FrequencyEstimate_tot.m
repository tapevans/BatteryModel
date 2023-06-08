function [re_est,re_est_var,im_est,im_est_var,nomsys]=FrequencyEstimate_tot(x,np,nz,P,w)
%FREQUENCYESTIMATE_TOT
%
%[re_est,re_est_var,im_est,im_est_var]=FrequencyEstimate_tot(x,np,nz,P,w)
%
% inputs: 
%       x - optimization vector from XXX
%       np - number of poles
%       nz - number of zeros
%       P - covariance of x
%       w - frequencies at which to estimate the impedance
%
% outputs: 
%       re_est - estimate of real impedance 
%       re_est_var - estimate variance
%       im_est - estimate of imaginary impedance
%       im_est_var - estimate variance



% nominal sysm
%
k = x(1);
p = x(1+(1:np));
z = x(1+np+(1:nz));
nomsys = zpk(z,p,k);
[renom,imnom] = nyquist(nomsys,w);
renom=squeeze(renom);
imnom=squeeze(imnom);


dFreqre_dparam=[];
dFreqim_dparam=[];
for ell=1:length(x)
    del = zeros(size(x));
    del(ell) = max(abs(x(ell))*1e-6,1e-9);
    dx = x-del;
    dk = dx(1);
    dp = dx(1+(1:np));
    dz = dx(1+np+(1:nz));
    delsys = zpk(dz,dp,dk);
    [redel,imdel] = nyquist(delsys,w);
    redel=squeeze(redel);
    imdel=squeeze(imdel);
    dFreqre_dparam(:,ell) = (redel-renom)/del(ell);
    dFreqim_dparam(:,ell) = (imdel-imnom)/del(ell);
end;


re_est = renom;
im_est = imnom;
%
% get estimated variance at each frequency
%
re_est_var = diag(dFreqre_dparam*P*dFreqre_dparam');
im_est_var = diag(dFreqim_dparam*P*dFreqim_dparam');
   