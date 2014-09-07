entier='int32';
reel='float32';
texte='char';
precision='ab';

type='bl';
byte_in='big_endian';
byte_out='little_endian';
fi=fopen('crsa1','r',type(1));
a=fread(fi,2,entier);
fclose(fi)
if abs(a(1))>1000
  type='lb';
  byte=byte_in;
  byte_in=byte_out;
  byte_out=byte;
  clear byte
end

dir_out=['dir.' byte_out];
eval(['!mkdir ' dir_out]);
for ific=1:8
  sfic=int2str(ific);
  for iprecision=1:2
    file_in=['crs' precision(iprecision) sfic];
    file_out=[dir_out '/' file_in];
    fin=fopen(file_in,'r',type(1));
    fout=fopen(file_out,'w',type(2));
    disp(['fichier ' file_in ' traite']);

%   nen,nspec,nbrexc,nbrionst
    a=fread(fin,7,entier);
    fwrite(fout,a,entier);
    disp(a(2:5)')
    nen=a(1);
    nspec=a(2);
    nbrexc=a(3);
    nbrionst=a(4);

%   (centE(n),n=1,nen)
    len=a(end)/4;
    a=fread(fin,len,reel);
    fwrite(fout,a,reel);
    em=a(1);
    a=fread(fin,2,entier);
    fwrite(fout,a,entier);

%   (botE(n),n=1,nen)
    len=a(end)/4;
    a=fread(fin,len,reel);
    eb=a(1);
    fwrite(fout,a,reel);
    a=fread(fin,2,entier);
    fwrite(fout,a,entier);

%   (engdd(n),n=1,nen)
    len=a(end)/4;
    a=fread(fin,len,reel);
    de=a(1);
    fwrite(fout,a,reel);
    a=fread(fin,2,entier);
    fwrite(fout,a,entier);
    disp([eb em de])

%   (nexcst(isp),isp=1,nspec),(nionst(isp),isp=1,nspec)
    len=a(end)/4+2;
    a=fread(fin,len,entier);
    fwrite(fout,a,entier);
%   ((title(js,isp),isp=1,nspec),js=1,nbrexc)
    len=a(end);
    a=fread(fin,len,texte);
    fwrite(fout,a,texte);
    a=fread(fin,2,entier);
    fwrite(fout,a,entier);

%   (((ethres(isp,js,jp),isp=1,nspec),js=1,nbrexc),jp=1,nbrionst)
    len=a(end)/4;
    a=fread(fin,len,reel);
    fwrite(fout,a,reel);
    a=fread(fin,2,entier);
    fwrite(fout,a,entier);

%   ((bratio(jp,isp),jp=1,nbrionst),isp=1,nspec)
    len=a(end)/4;
    a=fread(fin,len,reel);
    fwrite(fout,a,reel);
    a=fread(fin,2,entier);
    fwrite(fout,a,entier);

%   ((cel(j,n),j=1,nspec),n=1,nen)
    len=a(end)/4;
    a=fread(fin,len,reel);
    fwrite(fout,a,reel);
    a=fread(fin,2,entier);
    fwrite(fout,a,entier);

%   ((cin(j,n),j=1,nspec),n=1,nen)
    len=a(end)/4;
    a=fread(fin,len,reel);
    fwrite(fout,a,reel);
    a=fread(fin,2,entier);
    fwrite(fout,a,entier);

%   (((cinex(j,js,n),j=1,nspec),js=1,nbrexc),n=1,nen)
    len=a(end)/4;
    a=fread(fin,len,reel);
    fwrite(fout,a,reel);
    a=fread(fin,2,entier);
    fwrite(fout,a,entier);

%   (((cinex(j,js,n),j=1,nspec),js=1,nbrexc),n=1,nen)
    len=a(end)/4;
    a=fread(fin,len,reel);
    fwrite(fout,a,reel);

    a=fread(fin,1,entier);
    fwrite(fout,a,entier);

    fclose(fin);
    fclose(fout);

    file_in=['rdt' precision(iprecision) sfic];
    file_out=[dir_out '/' file_in];
    fin=fopen(file_in,'r',type(1));
    fout=fopen(file_out,'w',type(2));
    disp(['fichier ' file_in ' traite']);

%   line
    a=fread(fin,1,entier);
    fwrite(fout,a,entier);
    len=a(end);
    a=fread(fin,len,texte);
    fwrite(fout,a,texte);

%   nen,nspec
    a=fread(fin,4,entier);
    fwrite(fout,a,entier);
    disp(a(3:4)')
    nen=a(3);

    for i=1:nen

%   ((omdeg(n,nfix,isp),n=nen,nend,-1),isp=1,nspec)
      a=fread(fin,2,entier);
      fwrite(fout,a,entier);
      len=a(end)/4;
      a=fread(fin,len,reel);
      fwrite(fout,a,reel);
      a=fread(fin,2,entier);
      fwrite(fout,a,entier);

%   ((omdeg(nfix,n,isp),n=nen,nend,-1),isp=1,nspec)
      len=a(end)/4;
      a=fread(fin,len,reel);
      fwrite(fout,a,reel);

    end


    a=fread(fin,1,entier);
    fwrite(fout,a,entier);

    fclose(fin);
    fclose(fout);

  end
end

