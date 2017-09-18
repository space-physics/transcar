a=dir('*_*_*.mat');
filelist=[];
for i=1:length(a)
 filelist=[filelist;a(i).name];
end
nbf=size(filelist);nbf=nbf(1);

data_path=pwd;
data_type='tra';

tube_time=[];
tubelon_time=[];
tubelat_time=[];
tubetmag_time=[];
tubelonmag_time=[];
tubelatmag_time=[];
tubeEnord_time=[];
tubeEest_time=[];

save tube tube_time tubelon_time tubelat_time ...
     tubetmag_time  tubelonmag_time tubelatmag_time ...
     tubeEnord_time tubeEest_time
for ifile=1:nbf
  file=filelist(ifile,:)
  giveno
  close all
  load tube
  ntime=length(time);
  tube_time(1:ntime,ifile)=time(:);
  tubelon_time(1:ntime,ifile)=lon_time(:);
  tubelat_time(1:ntime,ifile)=lat_time(:);
  tubetmag_time(1:ntime,ifile)=tmag_time(:);
  tubelonmag_time(1:ntime,ifile)=lonmag_time(:);
  tubelatmag_time(1:ntime,ifile)=latmag_time(:);
  tubeEnord_time(1:ntime,ifile)=Enord_time(:);
  tubeEest_time(1:ntime,ifile)=Eest_time(:);

  save tube tube_time tubelon_time tubelat_time ...
     tubetmag_time  tubelonmag_time tubelatmag_time ...
     tubeEnord_time tubeEest_time

end

i=find(~tube_time(end,:));
tube_time(end,i)=tube_time(end-1,i);
tube_time(end,i)=tube_time(end-1,i);
tubelon_time(end,i)=tubelon_time(end-1,i);
tubelat_time(end,i)=tubelat_time(end-1,i);
tubetmag_time(end,i)=tubetmag_time(end-1,i);
tubelonmag_time(end,i)=tubelonmag_time(end-1,i);
tubelatmag_time(end,i)=tubelatmag_time(end-1,i);
tubeEnord_time(end,i)=tubeEnord_time(end-1,i);
tubeEest_time(end,i)=tubeEest_time(end-1,i);

save tube tube_time tubelon_time tubelat_time ...
     tubetmag_time  tubelonmag_time tubelatmag_time ...
     tubeEnord_time tubeEest_time

clear all
load tube
