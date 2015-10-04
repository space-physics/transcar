function big2little(filein,fileout);

fin=fopen(filein,'r','b');

if fin>0
  fout=fopen(fileout,'w','l');
  a=fread(fin,'int32');
  fwrite(fout,a,'int32');
  fclose(fin);
  fclose(fout);
else
  disp(['erreur d''ouverture du fichier' filein]')
end
