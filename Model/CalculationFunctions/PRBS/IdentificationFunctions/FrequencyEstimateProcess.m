function [re_est,re_est_var,im_est,im_est_var,re,im,var_re,var_im]=FrequencyEstimateProcess(sys,wtot,Update_Freq,PRBSLength)

%
% Initialize variables
%
re=zeros(length(wtot),length(sys));
im=zeros(length(wtot),length(sys));
sdre=zeros(length(wtot),length(sys));
sdim=zeros(length(wtot),length(sys));

for k=1:length(sys),
    %
    % Find which frequencies of interest are inside the proper region for
    % this model
    % 
%    windex = find(wtot > pi/(dat{k}.ts)/1000 & wtot < pi/(dat{k}.ts)/5);
%    windex = find(wtot > pi/(dat{k}.ts)/(length(dat{k}.OutputData)) & wtot < pi/(dat{k}.ts));
    windex = find(wtot > 2*pi*Update_Freq{k}/PRBSLength & wtot < pi*Update_Freq{k});
    w=wtot(windex);
    %
    % numerical differentiation from parameters to frequency response
    %
    p=getpvec(sys{k})';
    P=getcov(sys{k});
    N = length(p);
    n = floor(N/2);
%%    nomsys = d2c(tf([p(1:N/2) 0],[1 p(N/2+1:end)],sys{k}.ts));
%    nomsys = d2c(tf([p(1:n+1)],[1 p(n+2:end)],sys{k}.ts));
    nomsys = d2c(sys{k});
    delcsys = sys{k};
    [renom,imnom] = nyquist(nomsys,w);
    renom=squeeze(renom);
    imnom=squeeze(imnom);
    dFreqre_dparam=[];
    dFreqim_dparam=[];
    for ell=1:length(p),
        del = zeros(1,length(p));
        del(ell) = max(abs(p(ell))*1e-6,1e-9);
        pdel = p + del;
 %       delsys = d2c(tf([pdel(1:N/2) 0],[1 pdel(N/2+1:end)],sys{k}.ts));
        delcsys = setpvec(delcsys,pdel);
        delsys = d2c(delcsys);
        [redel,imdel] = nyquist(delsys,w);
        redel=squeeze(redel);
        imdel=squeeze(imdel);
        dFreqre_dparam(:,ell) = (redel-renom)/del(ell);
        dFreqim_dparam(:,ell) = (imdel-imnom)/del(ell);
    end;
    
    %
    % get estimated standard deviation at each frequency
    %
    var_re(:,k) = inf*ones(size(wtot));
    var_im(:,k) = inf*ones(size(wtot));
    var_re(windex,k) = diag(dFreqre_dparam*P*dFreqre_dparam');
    var_im(windex,k) = diag(dFreqim_dparam*P*dFreqim_dparam');
   
    
    re(:,k) = zeros(size(wtot));
    im(:,k) = zeros(size(wtot));
    re(windex,k)=renom;
    im(windex,k)=imnom;
end;

if length(sys)>1,
 % Weighted sum estimates
 re_weight = (1./var_re);
 re_est = re.*(re_weight);
 re_est = sum(re_est')./sum(re_weight');
 re_est_var = (1./sum(re_weight'));
 im_weight = (1./var_im);
 im_est = im.*(im_weight);
 im_est = sum(im_est')./sum(im_weight');
 im_est_var = (1./sum(im_weight'));
else
 re_est = re;
 re_est_var = var_re;
 im_est = im;
 im_est_var = var_im;
end;