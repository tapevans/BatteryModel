%
% Figure with model fit
%
figure('Position',[500 500 700 500])
ax1=axes('Position',[.25 .11 .65 .8]);
Ts_all={};
for i=1:length(smoothdat)
data.sys{i} = G;
Ts_all{i} = smoothdat{i}.Ts;
end;
data.smoothdat = smoothdat;
title(ax1,'Compare')
data.ax1=ax1;
aa=uicontrol('Style','popup','String',cellfun(@num2str,Ts_all, 'UniformOutput',false),'Position', [20 340 100 50],'UserData',data,'Callback',@EISresultsdisplaycompare_callback);
callbackCell = get(aa,'Callback');
feval(callbackCell,aa)

figure
hold on
uncertainty = max(re_est_var,im_est_var);
index=find(3*sqrt(uncertainty)<.001);

errorbar(1000*re_est(index),-1000*im_est(index),3*1000*sqrt(im_est_var(index)),3*1000*sqrt(im_est_var(index)),3*1000*sqrt(re_est_var(index)),3*1000*sqrt(re_est_var(index)),'*','color',[.8 .8 .8])
h=plot(1000*re_est(index),-1000*im_est(index),'b*');
    hold off
xlabel('Real Impedance (m \Omega)')
ylabel('-Imaginary Impedance (m \Omega)')
a = axis;
[dummy,ind]=min(wtot(index));
try
  hh=annotate(h,1000*re_est(ind),['f=',num2str(wtot(ind)/(2*pi)),' Hz'],'ur');
catch
end

[dummy,ind]=max(wtot(index));
try
  hh=annotate(h,1000*re_est(ind),['f=',num2str(wtot(ind)/(2*pi)),' Hz'],'lr');
catch
end;

ind=min(find(wtot<= sqrt(max(wtot)*min(wtot))));
try
  annotate(h,1000*re_est(ind),['f=',num2str(wtot(ind)/(2*pi)),' Hz'],'ul');
catch
end;
%title(savefilename)
set(gca,'fontsize',14)
a(2)=a(2)+1;
axis(a)
axis equal