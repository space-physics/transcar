fic=input('Nom du fichier : ','s');
eval(['load ' fic]);
k=findstr(fic,'.');
if ~isempty(k)
   fic=fic(1:(k-1));
end
eval(['trace=' fic ';']);
itube=trace(:,1);
t=trace(:,2)/3600;
tm=trace(:,3);
lm=trace(:,4);
lo=trace(:,5);
la=trace(:,6);
pot=trace(:,7);
ipos=trace(:,end);
clear trace
i=find(diff(itube)>0);
i=[0;i(:);length(itube)];
longueur=0;
for j=1:length(i)-1
   eval(['tps' int2str(j) '=t(i(j)+1:i(j+1));']);
   eval(['tm' int2str(j) '=tm(i(j)+1:i(j+1));']);
   eval(['lm' int2str(j) '=lm(i(j)+1:i(j+1));']);
   eval(['lo' int2str(j) '=lo(i(j)+1:i(j+1));']);
   eval(['la' int2str(j) '=la(i(j)+1:i(j+1));']);
   eval(['pot' int2str(j) '=pot(i(j)+1:i(j+1));']);
   eval(['he' int2str(j) '=fix(tps' int2str(j) ');']);
   eval(['mi' int2str(j) '=60*(tps' int2str(j) '-he' int2str(j) ');']);
   eval(['se' int2str(j) '=3600*(tps' int2str(j) '-he' int2str(j) '-mi' int2str(j) '/60);']);
   eval(['stl' int2str(j) '=rem(tps' int2str(j) '+lo' int2str(j) '/15+24,24);']);
   eval(['[Chi' int2str(j) ',zmin' int2str(j) '] = kiangle (lo' int2str(j) ',la' int2str(j) ',1993,2,16,he' int2str(j) ',mi' int2str(j) ',se' int2str(j) ');']);
   eval(['[x' int2str(j) ',y' int2str(j) ']=loc2cart(tm' int2str(j) ',lm' int2str(j) ');']);
   eval(['t_st(j)=tm' int2str(j) '(end);']);
   eval(['t_en(j)=tps' int2str(j) '(1);']);
   longueur=max(longueur,i(j+1)-i(j));
%   eval(['plot3(' int2str(j) '*ones(size(tps' int2str(j) ')),tps' int2str(j) ',Chi' int2str(j) ')']);
   
%   plot(x,y)
%   hold on
end
x=ones([longueur,length(i)-1])*NaN;
y=ones([longueur,length(i)-1])*NaN;
Chi=ones([longueur,length(i)-1])*NaN;
tps=ones([longueur,length(i)-1])*NaN;
stl=ones([longueur,length(i)-1])*NaN;
tm=ones([longueur,length(i)-1])*NaN;
lo=ones([longueur,length(i)-1])*NaN;
la=ones([longueur,length(i)-1])*NaN;
lm=ones([longueur,length(i)-1])*NaN;
pot=ones([longueur,length(i)-1])*NaN;

for j=1:length(i)-1
  eval(['x(1:length(x' int2str(j) '),' int2str(j) ')=x' int2str(j) ';']);
  eval(['y(1:length(y' int2str(j) '),' int2str(j) ')=y' int2str(j) ';']);
  eval(['Chi(1:length(Chi' int2str(j) '),' int2str(j) ')=Chi' int2str(j) ';']);
  eval(['tps(1:length(tps' int2str(j) '),' int2str(j) ')=tps' int2str(j) ';']);
  eval(['stl(1:length(stl' int2str(j) '),' int2str(j) ')=stl' int2str(j) ';']);
  eval(['tm(1:length(tm' int2str(j) '),' int2str(j) ')=tm' int2str(j) ';']);
  eval(['lo(1:length(lo' int2str(j) '),' int2str(j) ')=lo' int2str(j) ';']);
  eval(['la(1:length(la' int2str(j) '),' int2str(j) ')=la' int2str(j) ';']);
  eval(['lm(1:length(lm' int2str(j) '),' int2str(j) ')=lm' int2str(j) ';']);
  eval(['pot(1:length(pot' int2str(j) '),' int2str(j) ')=pot' int2str(j) ';']);
end
lo=rem(lo+360,360);
