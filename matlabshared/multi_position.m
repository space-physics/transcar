clear per pos position_time ipos lat_ref lon_ref indi
per.az=NaN;
per.el=NaN;
per.lon=NaN;
per.lat=NaN;
per.alt=NaN;
len=1;

display('multi_position: create variable "multi_pos"')



if ~exist('multi_pos','var')
   len=1;
   nb_pos=1;
   per.az=az;
   per.el=el;
   per.alt=z(1);
   per.lon=longeo;
   per.lat=latgeo;
   pos(1).index=1:length(time);
   position_time=ones(size(time));

else

   delta_lat=.3;
   delta_lon=1;
   delta_alt=1;
   lat_ref=lat_time(1,:); 
   lon_ref=lon_time(1,:); 
   alt_ref=z_time(45,:); 
   
   for i=1:length(lat_ref)
       j=find(abs([per.lat]-lat_ref(i))< delta_lat & abs([per.lon]-lon_ref(i))<=delta_lon);
       if isempty(j)
          per(len).lat=lat_ref(i);
          per(len).lon=lon_ref(i);
          per(len).az= az_time(i);
          per(len).el= el_time(i);
          per(len).alt=alt_ref(i);
          len=len+1;
       end
   end
   len=len-1;
   latit=[per.lat];
   longi=[per.lon];
   elev =[per.el];
   azim =[per.az];
   altit=[per.alt];
   [latit ilatsort]=sort(latit);
   altit=altit(ilatsort);
   longi=longi(ilatsort);
   elev =elev (ilatsort);
   azim =azim (ilatsort);
   for i=1:len
      per(i).lat=latit(i);
      per(i).lon=longi(i);
      per(i).az =azim (i);
      per(i).el =elev (i);
      per(i).alt=altit(i);
   end

   for i=len:-1:1
      pos(i).index=find(abs(lat_ref-per(i).lat)< delta_lat & abs(lon_ref-per(i).lon)<=delta_lon);
      indi(i)=length(pos(i).index);
   end
   nb_pos=len;
   position_time=ones(size(time))*NaN;
   for j=1:nb_pos
     position_time(pos(j).index)=j;
   end

   scan_num=1;
   scan=ones(size(time))*NaN;
   scan_deb=1;
   while isnan(position_time(scan_deb)==1)
      scan_deb=scan_deb+1;
   end
   scan(scan_deb)=scan_num;
   last_pos=position_time(scan_deb);
   last_time=time(scan_deb);
   % Dans ce qui suit, on suppose qu on va du nord vers le sud; 
   % comme le scan st ordonne en latitude croissantes, le sens "normal de scan est en 
   % positions decroissantes. Si c'est l'inverse comme pour cp3-e (je crois), alors 
   % il faut remplacer la ligne ci-dessous "if position_time(i)<=last_pos..." par:
   %                                       "if position_time(i)>=last_pos..."
   
   for i=scan_deb+1:length(time)
      if isnan(position_time(i))==0
         if position_time(i)<=last_pos & time(i)<last_time+0.5
            last_time=time(i);
            last_pos=position_time(i);
            scan(i)=scan_num;
         elseif position_time(i)>0
            last_time=time(i);
            last_pos=position_time(i);
            scan_num=scan_num+1;
            scan(i)=scan_num;
         else
            scan(i)=NaN;
         end
      end
   end
   if nb_pos > 1
      j=diff(position_time);
      i=find(j==0);
      position_time(i)=NaN;
      scan(i)=NaN;
   end
end
