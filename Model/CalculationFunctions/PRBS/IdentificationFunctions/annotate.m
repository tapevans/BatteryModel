function hh=annotate(h,pos,texttag,direction,angle,scale,yscale,fontsize,interpreter)
%ANNOTATE(h,pos,texttag,direction,angle,scale,yscale,fontsize,interpreter) add the an annotation to a plot
%
% h = handle of plot
% pos = x position for annotation
% texttag = text of annotation
% direction = 'lr' lower right (default)
%           = 'ur' upper right
%           = 'll' lower left
%           = 'ul' upper left
% angle = angle of line (default pi/6)
% scale = (default .05)
% yscale = (default 1)
% fontsize = fontsize (default 14)
% interpreter = 1 for latex (default 1)
if (nargin<4 | isempty(direction)), direction='lr'; end;
if (nargin<5 | isempty(angle)), angle=pi/6; end;
if (nargin<6 | isempty(scale)), scale=.05; end;
if (nargin<7 | isempty(yscale)), yscale=1; end;
if (nargin<8 | isempty(fontsize)), fontsize=14; end;
if (nargin<9 | isempty(interpreter)), interpreter=1; end;

xdata=get(h,'xdata');
[xdata,i]=sort(xdata);
ydata=get(h,'ydata');
ydata=ydata(i);
xdiff=diff(xdata);
ni=find(xdiff==0);
i=setdiff(1:length(xdata),ni);
xdata=xdata(i);
ydata=ydata(i);

aa=get(h,'parent');
ax=get(aa,'xlim');
ay=get(aa,'ylim');
switch direction;
    case 'ur'
      theta=angle;
      halign='Left';
      len=scale;
    case 'lr'
      theta=-angle;
      halign='Left';
      len=scale;
    case 'ul'
      theta=-angle;
      halign='Right';
      len=-scale;
    case 'll'
      theta=angle;
      halign='Right';
      len=-scale;
    otherwise
       error('Invalid direction');
end;

if pos>max(xdata) | pos<min(xdata), warning('Position of annotation not within max and min of plot data'); end;

if (strcmp(get(gca,'Xscale'),'linear')),
 dist=[len*cos(theta)*(ax(2)-ax(1)), yscale*len*sin(theta)*(ay(2)-ay(1))];
 p1=[pos interp1(xdata,ydata,pos)];
 p2=p1+dist; 
 p3=p2+[len*(ax(2)-ax(1)) 0];
elseif (strcmp(get(gca,'Xscale'),'log')),
 dist=[len*cos(theta)*(log10(ax(2))-log10(ax(1))), yscale*len*sin(theta)*(ay(2)-ay(1))];  
 p1=[pos interp1(xdata,ydata,pos)];
 p2=[10^(log10(p1(1))+dist(1)),p1(2)+dist(2)]; 
 p3=[10^(log10(p2(1))+len*(log10(ax(2))-log10(ax(1)))),p2(2)];
else
    error('Unsupported axis scaling')
end;
hh(1)=line([p1(1) p2(1) p3(1)],[p1(2) p2(2) p3(2)],'linewidth',1);
if interpreter==1,
  hh(2)=text(p3(1),p3(2),texttag,'HorizontalAlignment',halign,'VerticalAlignment','Middle','Fontsize',fontsize,'interpreter','latex');
else
  hh(2)=text(p3(1),p3(2),texttag,'HorizontalAlignment',halign,'VerticalAlignment','Middle','Fontsize',fontsize);
end;    
