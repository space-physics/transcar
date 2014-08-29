function create_precinput(tset,er,phimem,filelab)
%% MAKE SURE ENERGY GRID IS IN BINS WITH BOTH A BEGINNING AND ENDING ENERGY
[ne,nt]=size(phimem);
if length(er)==ne
   de=diff(er);
   er=[er(:);er(ne)+de(length(de))];
end
%% OPEN THE FILE AND WRITE TO IT
fid=fopen([filelab,'.dat'],'w');
for k=1:nt
    fprintf(fid,'%f\n',round(tset(k)*3600));
    data=[er(:),[phimem(:,k);-1.0]]';
    fprintf(fid,'%e %e\n',data);
end
fprintf(fid,'%f\n',round((tset(nt)+(tset(nt)-tset(nt-1)))*3600));
fprintf(fid,'-1.0 -1.0\n');
fclose(fid);
end