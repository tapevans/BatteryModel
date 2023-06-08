function EISresultsdisplayEIS_callback(source,event)
data = source.UserData;
for i=1:length(data.handles)
           index(i)=data.handles{i}.Value;
end;

re = data.re;
var_re = data.var_re;
im = data.im;
var_im = data.var_im;
ax1 = data.ax1;

 plot(ax1,re(~isinf(var_re(:,1)),1),-im(~isinf(var_re(:,1)),1),'^','color',[.9 .9 .9])
hold on
 for i=2:size(re,2),
     plot(ax1,re(~isinf(var_re(:,i)),i),-im(~isinf(var_re(:,i)),i),'^','color',[.9 .9 .9])
 end;


index=logical(index);

re=re(:,index);
im=im(:,index);
var_re=var_re(:,index);
var_im=var_im(:,index);


if size(re,2)>1,
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
re_est = re_est(~isinf(re_est_var));
re_est_var = re_est_var(~isinf(re_est_var));
im_est = im_est(~isinf(im_est_var));
im_est_var = im_est_var(~isinf(im_est_var));

plot(ax1,re_est,-im_est,'*')

%
%errorbar(ax1,re_est,-im_est,3*sqrt(im_est_var),3*sqrt(im_est_var),3*sqrt(re_est_var),3*sqrt(re_est_var),'*')
    hold off
axis equal