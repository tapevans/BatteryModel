%
% Figure with impedances
%
figure('Position',[500 10 700 900])
% figure
ax1=axes('Position',[.25 .5 .65 .45]);
ax2=axes('Position',[.25 .11 .65 .3]);
data.re = re;
data.im = im;
data.var_re = var_re;
data.var_im = var_im;
data.ax1=ax1;
data.ax2=ax2;
data.sys=sys;
data.wtot = wtot;
title(ax1,'EIS plot')
w=uicontrol('Style','edit','String','[.03 .06 0 .01]','Position', [20 540 80 30]);
data.w=w;
aa=uicontrol('Style','popup','String',cellfun(@num2str,Ts_all, 'UniformOutput',false),'Position', [20 640 100 50],'UserData',data,'Callback',@EISresultsdisplayfreq_callback);
callbackCell = get(aa,'Callback');
feval(callbackCell,aa)


%
% Figure with model fit
%
figure('Position',[500 10 700 500])
ax1=axes('Position',[.25 .11 .65 .8]);
data.sys = sys;
data.smoothdat = smoothdat;
title(ax1,'Compare')
data.ax1=ax1;
aa=uicontrol('Style','popup','String',cellfun(@num2str,Ts_all, 'UniformOutput',false),'Position', [20 340 100 50],'UserData',data,'Callback',@EISresultsdisplaycompare_callback);
callbackCell = get(aa,'Callback');
feval(callbackCell,aa)

figure('Position',[500 10 700 500])
ax1=axes('Position',[.25 .11 .65 .8]);
data.re = re;
data.var_re = var_re;
data.im = im;
data.var_im = var_im;
title(ax1,'EIS')
data.ax1=ax1;
aa={};
for i=1:length(Ts_all),
    aa{i}=uicontrol('Style','checkbox','Value',1,'String',num2str(Ts_all{i}),'Position',[20 340-50*i 100 50],'Callback',@EISresultsdisplayEIS_callback);
end;
data.handles = aa;
for i=1:length(Ts_all),
    set(aa{i},'UserData',data)
end;
callbackCell = get(aa{1},'Callback');
feval(callbackCell,aa{1})