

offset=(rec_num-1)*size_record*int_type;
status=fseek(fid,offset,'bof');
a=fread(fid,size_record,type);

b=reshape(a,ncol,nli);

head=b(1:ncol,1:2);
data=b(1:ncol,3:nli)';

an =head(3);
mo =head(4);
jo =head(5);
he =head(6);
mi =head(7);
se =head(8);
intpas=head(9);
%temps=datenum(an,mo,jo,he,mi,se)+intpas/86400/2;
%[an,mo,jo,he,mi,se]=datevec(temps);

longeo   =head(10);
latgeo   =head(11);
lonmag   =head(12);
latmag   =head(13);
tmag     =head(14);
f107     =head(15);
f107a    =head(16);
ap       =head(17);
kp       =head(18);
dtn      =head(19);
dun      =head(20);
cofo     =head(21);
cofh     =head(22);
cofn     =head(23);
chi      =head(24);
Fe       =head(25);
Ee       =head(26);
Fi       =head(27);
Ei       =head(28);
Bmag     =head(29);
dipangle =head(30);
Enord    =head(31);
Eest     =head(32);
vperpnord=head(33);
vperpest =head(34);
vhorizon =head(35);
vpara    =head(36);
approx   =head(37);
ddp      =head(38);
Jtop     =head(39);

ipos_z=1;
ipos_n1=2;
ipos_n2=3;
ipos_n3=4;
ipos_n4=5;
ipos_n5=6;
ipos_n6=7;
ipos_v1=8;
ipos_v2=9;
ipos_v3=10;
ipos_vm=11;
ipos_ve=12;
if approx>=13
  ipos_t1p=13;
  ipos_t1t=14;
  ipos_t2p=15;
  ipos_t2t=16;
  ipos_t3p=17;
  ipos_t3t=18;
  ipos_tmp=19;
  ipos_tmt=20;
  ipos_tep=21;
  ipos_tet=22;
  ipos_q1=23;
  ipos_q2=24;
  ipos_q3=25;
  ipos_qe=26;
  ipos_nno=27;
  ipos_uno=28;
  ipos_po=29;
  ipos_ph=30;
  ipos_pn=31;
  ipos_pn2=32;
  ipos_po2=33;
  ipos_heat=34;
  ipos_no=35;
  ipos_nh=36;
  ipos_nn=37;
  ipos_nn2=38;
  ipos_no2=39;
  ipos_tn=40;
  ipos_un=41;
  ipos_vn=42;
  ipos_wn=43;
  ipos_nes=44;
  ipos_jes=45;
  ipos_tes=46;
  ipos_qes=47;
else
  ipos_t1p=13;
  ipos_t1t=13;
  ipos_t2p=14;
  ipos_t2t=14;
  ipos_t3p=15;
  ipos_t3t=15;
  ipos_tmp=16;
  ipos_tmt=16;
  ipos_tep=17;
  ipos_tet=17;
  ipos_q1=18;
  ipos_q2=19;
  ipos_q3=20;
  ipos_qe=21;
  ipos_nno=22;
  ipos_uno=23;
  ipos_po=24;
  ipos_ph=25;
  ipos_pn=26;
  ipos_pn2=27;
  ipos_po2=28;
  ipos_heat=29;
  ipos_no=30;
  ipos_nh=31;
  ipos_nn=32;
  ipos_nn2=33;
  ipos_no2=34;
  ipos_tn=35;
  ipos_un=36;
  ipos_vn=37;
  ipos_wn=38;
  if ncol>38
    ipos_nes=39;
    ipos_jes=40;
    ipos_tes=41;
    ipos_qes=42;
  else
    ipos_nes=ipos_z;
    ipos_jes=ipos_z;
    ipos_tes=ipos_z;
    ipos_qes=ipos_z;
  end
end


z =data(:,ipos_z);
n1=data(:,ipos_n1);
n2=data(:,ipos_n2);
n3=data(:,ipos_n3);
n4=data(:,ipos_n4);
n5=data(:,ipos_n5);
n6=data(:,ipos_n6);
v1=data(:,ipos_v1);
v2=data(:,ipos_v2);
v3=data(:,ipos_v3);
vm=data(:,ipos_vm);
ve=data(:,ipos_ve);
t1p=data(:,ipos_t1p);
t1t=data(:,ipos_t1t);
t1=(t1p+2*t1t)/3;
t2p=data(:,ipos_t2p);
t2t=data(:,ipos_t2t);
t2=(t2p+2*t2t)/3;
t3p=data(:,ipos_t3p);
t3t=data(:,ipos_t3t);
t3=(t3p+2*t3t)/3;
tmp=data(:,ipos_tmp);
tmt=data(:,ipos_tmt);
tm=(tmp+2*tmt)/3;
tep=data(:,ipos_tep);
tet=data(:,ipos_tet);
te=(tep+2*tet)/3;
q1=data(:,ipos_q1);
q2=data(:,ipos_q2);
q3=data(:,ipos_q3);
qe=data(:,ipos_qe);
Nno=data(:,ipos_nno);
Vno=data(:,ipos_uno);

nm=n4+n5+n6;
ne=n1+n2+n3+nm;
vion=(n1.*v1+n2.*v2+n3.*v3+nm.*vm)./ne;

Po=data(:,ipos_po);
Ph=data(:,ipos_ph);
Pn=data(:,ipos_pn);
Pn2=data(:,ipos_pn2);
Po2=data(:,ipos_po2);
Heat=data(:,ipos_heat);
No=data(:,ipos_no);
Nh=data(:,ipos_nh);
Nn=data(:,ipos_nn);
Nn2=data(:,ipos_nn2);
No2=data(:,ipos_no2);
Tn=data(:,ipos_tn);
Un=data(:,ipos_un);
Vn=data(:,ipos_vn);
Wn=data(:,ipos_wn);

if ncol>38,
  nes=data(:,ipos_nes);
  jes=data(:,ipos_jes);
  tes=data(:,ipos_tes);
  qes=data(:,ipos_qes);
end

az=head(ncol,2);
el=dipangle;

%Added for nuetral Hot O
if ncol>47
   NOHot=data(:,48);
   TnOHot=data(:,49);
   
   if ncol>49
     n7=data(:,50);
     v7=data(:,51);
     t7p=data(:,52);
     t7t=data(:,53);
     t7=(t7p+2*t7t)/3;
     q7=data(:,54);
     POHot=data(:,55);
     
     ne=ne+n7;
     
     if ncol>55
        nOHot=data(:,56);
        vOHot=data(:,57);
        tOHotp=data(:,58);
        tOHott=data(:,59);
        tOHot=(tOHotp+2*tOHott)/3;
        qOHot=data(:,60);
        
        if ncol>60
           Po1d=data(:,61);
           no1d=data(:,62);
           vo1d=data(:,63);
        end
     end
   end
end
%-MZ

