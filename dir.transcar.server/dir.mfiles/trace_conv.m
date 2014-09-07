%temps=0:.5:24;
%long=temps*15-90;
%lat=60:90;
%[ax,ay,az]=pol2geo(long,lat);
%[flux,fluxenerg]=hardelec(lat,temps,3);
%flux=10.^flux;
%fluxenerg=10.^fluxenerg;
%energ=fluxenerg./flux;
%fluxenerg=fluxenerg*4*pi*1.6e-9;
%flux=flux*4*pi*1.e4;

[x,y]=loc2cart(tubetmag_time,tubelatmag_time);

[is,it]=size(tube_time);
nb_pos=10;
nb_temps=it/nb_pos;

ts=ceil(3600*max(tube_time(1,:)));
tf=floor(min(tube_time(end,:))*3600);
ts=min(tube_time(1,:))*3600;
tf=max(tube_time(end,:))*3600;
T=ts:60:tf;
T=T/3600;
clear X Y Xscan Yscan

for i=1:it
  for j=1:length(T)
    [a,ii]=min(abs(tube_time(:,i)-T(j)));
    X(j,i)=x(ii,i);
    Y(j,i)=y(ii,i);
  end
end
i_scan=0:(nb_temps-1);
j_scan=5+i_scan;
ind_scan=j_scan+i_scan*nb_pos;
k=find(j_scan<=nb_pos);
ind_scan=ind_scan(k);
i_scan=i_scan(k);
j_scan=j_scan(k);
it=length(ind_scan);
for i=1:it
  for j=1:length(T)
    [a,ii]=min(abs(tube_time(:,ind_scan(i))-T(j)));
    Xscan(j,i)=x(ii,ind_scan(i));
    Yscan(j,i)=y(ii,ind_scan(i));
  end
end

for i=1:length(T)
%  plot(x,y,x([1 end],:),y([1 end],:),'rp',X(i,:)',Y(i,:),'b.');
%  pcolor(ax,ay,log10(energ));colormap(jet);shading interp
%  hold on
  if (i>=32)
    plot(X(i,:)',Y(i,:),'b.',Xscan(i,:)',Yscan(i,:),'r*');
  else
    plot(X(i,:)',Y(i,:),'b.');
  end
  if (i>=32&i<=36)
    hold on
    h=plot(Xscan(i,:)',Yscan(i,:),'y');
    set(h,'linewidth',3);
  end
%  cadre([0 360],[70 90],'g')
   cadre([90 180],[70 90],'g');axis([-.4 .01 -.01 .4])
  set(gca,'visible','off');
  axis('square')
  idate=strrep(sprintf('%2d',i),' ','0');
  disp(idate)
  if i==length(T)
    hold on
%    plot([x(end,1:nb_pos:end);x(end,nb_pos:nb_pos:end)],[y(end,1:nb_pos:end);y(end,nb_pos:nb_pos:end)],'m')
    h=plot(x(end,ind_scan),y(end,ind_scan),'y');
    set(h,'linewidth',3);
  end
%  eval(['print -djpeg90 fig' idate]);
  pause(.1)
end