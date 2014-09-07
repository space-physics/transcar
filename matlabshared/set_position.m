clear('index')

  ipos=find(position_time==req_pos);
  
if nb_records>1
  dt=diff(time(ipos));
  dtmin=mean(dt);
  [dt,indt]=find(dt>1.5*dtmin);
  indt=[0;indt(:);length(ipos)];
  nin=length(indt);
  for i=1:nin-1
    index(i).data=ipos(indt(i)+1:indt(i+1));
  end
  

end
%%
if nb_pos>1
   if remote || strcmp(data_type(1:3),'tra')
      string_box=sprintf(' Lat=%2.1f Lon=%2.1f',...
                         per(req_pos).lat,per(req_pos).lon);
   else
      string_box=sprintf(' Az=%3d El=%3d',...
                         round(per(req_pos).az),round(per(req_pos).el));
   end
   titre=[file_expr string_box];
else
   titre=file_expr;                       
end                     

titre=strrep(titre,'_','\_');
%%
if remote
   sort_alt;
   z_sel=z_choice(length(z_choice));
else
   diff_z=diff(z_time(ipos(1)));
   sup_dfz=find(diff_z<0);
   if ~isempty(sup_dfz)
     ind_dfz=length(sup_dfz)+1;
   else
     ind_dfz=1;
   end
   sup_dfz(ind_dfz)=npt;
   inf_dfz(1)=1;
   inf_dfz(2:ind_dfz)=sup_dfz(1:(ind_dfz-1))+1;
end
