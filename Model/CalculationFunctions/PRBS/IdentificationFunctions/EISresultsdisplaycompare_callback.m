function EISresultsdisplayfreq_callback(source,event)
k=source.Value;
data = source.UserData;
smoothdat = data.smoothdat;
sys = data.sys;
% [YH, FIT, X0]= compare(sys{k},smoothdat{k})
compare(sys{k},smoothdat{k})