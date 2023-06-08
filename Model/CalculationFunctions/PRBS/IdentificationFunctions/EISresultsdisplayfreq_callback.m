function EISresultsdisplayfreq_callback(source,event)
k=source.Value;
data = source.UserData;
re= data.re;
im = data.im;
var_re = data.var_re;
var_im = data.var_im;
ax1 = data.ax1;
ax2 = data.ax2;
sys = data.sys;
sysc={};
for i=1:length(sys)
    sysc{i} = d2c(sys{i});
end;
w = data.w;
wtot = data.wtot;
errorbar(ax1,re(:,k),-im(:,k),3*sqrt(var_re(:,k)),3*sqrt(var_re(:,k)),3*sqrt(var_im(:,k)),3*sqrt(var_im(:,k)),'*')
axis(ax1,eval(w.String))  
xlabel(ax1,'Real')
  ylabel(ax1,'Imaginary')
  title(ax1,'EIS plot')
axes(ax2)
bodeplot(sysc{k},wtot(~isinf(var_re(:,k))),'*')
hold on
bodeplot(sysc{:},{min(wtot), max(wtot)})
hold off