if exist('data_path','var')
  filein=[data_path,'dir.output/',file];
else
  filein=file;
end    

file_expr=['TRANSCAR - ',filein];

int_type=4;
type=['float',int2str(8*int_type)];

TCoutEndianness='b';
[fid, message]=fopen(filein,'r',TCoutEndianness);
if fid <= 2 %FIXME: is this correct?
  disp(['Error in opening file: ' filein])
  disp(message)
end

a=fread(fid,2,type);
fclose(fid);
if a(2)<1 %FIXME knowing what the data should be, detect endianness (seems risky)
  TCoutEndianness='l';
end

switch TCoutEndianness
    case 'l', TCtxt = 'little';
    case 'b', TCtxt = 'big';
end

display([filein,' being read as ',TCtxt,'-endian'])

[fid, message]=fopen(filein,'r',TCoutEndianness);
if fid <= 2
  display(['Error in opening file: ' filein])
  display(message)
end


a=fread(fid,2,type);
npt=a(1);
nli=npt+2;
ncol=a(2);
size_record=nli*ncol;
fseek(fid,0,'eof');
nb_records=fix(ftell(fid)/size_record/int_type);

%fclose(fid); %will be used in data_tra 
%fopen(filein,'r',TCoutEndianness); %mistake?

time=zeros([1 nb_records]);
an_time=time;
mo_time=time;
jo_time=time;
he_time=time;
mi_time=time;
se_time=time;

az_time=zeros([1 nb_records]);
el_time=zeros([1 nb_records]);
lon_time=zeros([1 nb_records]);
lat_time=zeros([1 nb_records]);
Fe_time=zeros([1 nb_records]);
Ee_time=zeros([1 nb_records]);
Fi_time=zeros([1 nb_records]);
Ei_time=zeros([1 nb_records]);
Eest_time=zeros([1 nb_records]);
Enord_time=zeros([1 nb_records]);
vpara_time=zeros([1 nb_records]);
vperpest_time=zeros([1 nb_records]);
vperpnord_time=zeros([1 nb_records]);
ddp_time=zeros([1 nb_records]);
Jtop_time=zeros([1 nb_records]);

z_time=zeros([npt nb_records]);
ne_time=z_time;
v1_time=z_time;
v2_time=z_time;
v3_time=z_time;
vm_time=z_time;
ve_time=z_time;
t1p_time=z_time;
t1t_time=z_time;
t1_time=z_time;
t2p_time=z_time;
t2t_time=z_time;
t2_time=z_time;
t3p_time=z_time;
t3t_time=z_time;
t3_time=z_time;
tmp_time=z_time;
tmt_time=z_time;
tm_time=z_time;
tep_time=z_time;
tet_time=z_time;
te_time=z_time;
chi_time=time;
tmag_time=time;
lonmag_time=time;
latmag_time=time;
Bmag_time=time;
dipangle_time=time;
Po_time=z_time;
Po2_time=z_time;
Pn2_time=z_time;
prod_time=z_time;
chauf_time=z_time;
Nno_time=z_time;
Vno_time=z_time;
Tn_time=z_time;
Nh_time=z_time;
No_time=z_time;
No2_time=z_time;
Nn2_time=z_time;
n1_time=z_time;
n2_time=z_time;
n3_time=z_time;
n4_time=z_time;
n5_time=z_time;
n6_time=z_time;
nm_time=z_time;
qe_time=z_time;
q1_time=z_time;
q2_time=z_time;
q3_time=z_time;
Un_time=z_time;
Vn_time=z_time;
Wn_time=z_time;

if ncol>38,
  nes_time=z_time;
  jes_time=z_time;
  tes_time=z_time;
  qes_time=z_time;
end

%Added for neutral Hot O
if ncol>47
    NOHot_time=z_time;
    TnOHot_time=z_time;
    
    if ncol>49
      n7_time=z_time;
      v7_time=z_time;
      t7p_time=z_time;
      t7t_time=z_time;
      t7_time=z_time;
      q7_time=z_time;
      POHot_time=z_time;
      
      if ncol>55
          nOHot_time=z_time;
          vOHot_time=z_time;
          tOHotp_time=z_time;
          tOHott_time=z_time;
          tOHot_time=z_time;        
          qOHot_time=z_time;
          
          if ncol>60
             Po1d_time=z_time;
             no1d_time=z_time;
             vo1d_time=z_time;
          end
      end
    end

    Heat_time=z_time;
end
%-MZ

dhe=0;
he0=0;

for rec_num=1:nb_records

  if mod(rec_num,2000)==0
    display(['read_tra: Processing record ',int2str(rec_num),...
              '  ; remaining ',int2str(nb_records-rec_num),' records'])
  end
  
  data_tra

  if he<he0,
    dhe=dhe+24;
  end
      
  he0=he;

  time(rec_num)=he+dhe+mi/60+se/3600;
%  time(rec_num)=temps;
  an_time(rec_num)=an;
  mo_time(rec_num)=mo;
  jo_time(rec_num)=jo;
  he_time(rec_num)=he;
  mi_time(rec_num)=mi;
  se_time(rec_num)=se;

  z_time(:,rec_num)=z;
  az_time(rec_num)=az;
  el_time(rec_num)=el;
  lat_time(rec_num)=latgeo;
  lon_time(rec_num)=longeo;
  chi_time(rec_num)=chi;
  tmag_time(rec_num)=tmag;
  latmag_time(rec_num)=latmag;
  lonmag_time(rec_num)=lonmag;
  Bmag_time(rec_num)=Bmag;
  dipangle_time(rec_num)=dipangle;
  Fe_time(rec_num)=Fe;
  Ee_time(rec_num)=Ee;
  Fi_time(rec_num)=Fi;
  Ei_time(rec_num)=Ei;
  Eest_time(rec_num)=Eest;
  Enord_time(rec_num)=Enord;
  vpara_time(rec_num)=vpara;
  vperpest_time(rec_num)=vperpest;
  vperpnord_time(rec_num)=vperpnord;
  ddp_time(rec_num)=ddp;
  Jtop_time(rec_num)=Jtop;
  
  ne_time(:,rec_num)=ne;
  n1_time(:,rec_num)=n1;
  n2_time(:,rec_num)=n2;
  n3_time(:,rec_num)=n3;
  n4_time(:,rec_num)=n4;
  n5_time(:,rec_num)=n5;
  n6_time(:,rec_num)=n6;
  nm_time(:,rec_num)=nm;
  if approx>=13
    tep_time(:,rec_num)=tep;
    tet_time(:,rec_num)=tet;
    t1p_time(:,rec_num)=t1p;
    t1t_time(:,rec_num)=t1t;
    t2p_time(:,rec_num)=t2p;
    t2t_time(:,rec_num)=t2t;
    t3p_time(:,rec_num)=t3p;
    t3t_time(:,rec_num)=t3t;
    tmp_time(:,rec_num)=tmp;
    tmt_time(:,rec_num)=tmt;
  else
    te_time(:,rec_num)=tep;
    t1_time(:,rec_num)=t1p;
    t2_time(:,rec_num)=t2p;
    t3_time(:,rec_num)=t3p;
    tm_time(:,rec_num)=tmp;
  end
  v1_time(:,rec_num)=v1;
  v2_time(:,rec_num)=v2;
  v3_time(:,rec_num)=v3;
  vm_time(:,rec_num)=vm;
  ve_time(:,rec_num)=ve;
  qe_time(:,rec_num)=qe;
  q1_time(:,rec_num)=q1;
  q2_time(:,rec_num)=q2;
  q3_time(:,rec_num)=q3;

  Po_time(:,rec_num)=Po;
  Po2_time(:,rec_num)=Po2;
  Pn2_time(:,rec_num)=Pn2;
  prod_time(:,rec_num)=Po+Po2+Pn2;
  chauf_time(:,rec_num)=Heat;
  Nno_time(:,rec_num)=Nno;
  Vno_time(:,rec_num)=Vno;
  Tn_time(:,rec_num)=Tn;
  Nh_time(:,rec_num)=Nh;
  No_time(:,rec_num)=No;
  No2_time(:,rec_num)=No2;
  Nn2_time(:,rec_num)=Nn2;
  Vn_time(:,rec_num)=Vn;
  Un_time(:,rec_num)=Un;
  Wn_time(:,rec_num)=Wn;
  
  %Added for neutral hot O
  if ncol>47
    NOHot_time(:,rec_num)=NOHot;
    TnOHot_time(:,rec_num)=TnOHot;
    
    if ncol>49
      n7_time(:,rec_num)=n7;
      v7_time(:,rec_num)=v7;
      t7p_time(:,rec_num)=t7p;
      t7t_time(:,rec_num)=t7t;
      t7_time(:,rec_num)=t7;
      q7_time(:,rec_num)=q7;
      POHot_time(:,rec_num)=POHot;
            
      if ncol>55
          nOHot_time(:,rec_num)=nOHot;
          vOHot_time(:,rec_num)=vOHot;
          tOHotp_time(:,rec_num)=tOHotp;
          tOHott_time(:,rec_num)=tOHott;
          tOHot_time(:,rec_num)=tOHot;
          qOHot_time(:,rec_num)=qOHot;
          
          if ncol>60
             Po1d_time(:,rec_num)=Po1d;
             no1d_time(:,rec_num)=no1d;
             vo1d_time(:,rec_num)=vo1d;
          end
      end
    end
   
     Heat_time(:,rec_num)=Heat;
  end
  %-MZ  
  
  if exist('nes_time'),
    nes_time(:,rec_num)=nes;
    jes_time(:,rec_num)=jes;
    tes_time(:,rec_num)=tes;
    qes_time(:,rec_num)=qes;
  end

end

j_time=(n1_time.*v1_time...
       +n2_time.*v2_time...
       +n3_time.*v3_time...
       +nm_time.*vm_time...
       -ne_time.*ve_time);

 if ncol>49
    j_time=j_time+n7_time.*v7_time;
 end
       
if exist('jes_time'),
%  j_time=j_time-jes_time;
end

j_time=1.6e-13*j_time;


if approx>=13
  t1_time=(t1p_time+2*t1t_time)/3;
  t2_time=(t2p_time+2*t2t_time)/3;
  t3_time=(t3p_time+2*t3t_time)/3;
  tm_time=(tmp_time+2*tmt_time)/3;
  te_time=(tep_time+2*tet_time)/3;
else
  t1p_time=t1_time;
  t1t_time=t1_time;
  t2p_time=t2_time;
  t2t_time=t2_time;
  t3p_time=t3_time;
  t3t_time=t3_time;
  tmp_time=tm_time;
  tmt_time=tm_time;
  tep_time=te_time;
  tet_time=te_time;
end
%% run multi_position
multi_position
