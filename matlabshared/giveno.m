% This graphic routine is able to plot EISCAT data save with
%
% 		    NCAR   format
% 		or GUISDAP format
% 		or LEJEUNE format
%
% or to plot results from the TRANSCAR ionospheric model
%
% It requires the definition of the input format in the string data_type :
%
%		data_type = 'ncar'
%		data_type = 'guisdap'
%		data_type = 'lejeune'
%		data_type = 'hygraph'
% or
%		data_type = 'transcar'
%
% Moreover the following variable strings should be defined
%
%		data_path = path to reach the file to be read
%
%		file      = name of the file to be read
%
% data_path may be a null string if file contains the full path. For
% GUISDAP format, file does not need to be defined.

deg2rad=pi/180;
%rad2deg=1./deg2rad;

%close all
%fclose('all');
clear('time','T','TL','w_GIVEME','w_alt_time','w_alt','w_myplot','w_time',...
      'dne','dn1','dnm','dte','dt1','dv1','dne_time','dEnord dEest',...
      'dn1_time','dnm_time','dte_time','dt1_time','dv1_time',...
      'dEnord_time','dEest_time','dzne','dzte','dzt1','dzv1','*_time')
format compact
req_time=0;
remote=0;
nbdigit=6;

time=[];
an_time=[];
mo_time=[];
jo_time=[];
he_time=[];
mi_time=[];
se_time=[];

z_time=[];
ne_time=[];
n1_time=[];
te_time=[];
t1_time=[];
v1_time=[];

az_time=[];
el_time=[];
az_sel=[];
el_sel=[];
PowerProfile=[];
Cp6_exp=[];
Cp4_exp=[];
Site_Number=[];

%ht=fclose('all');
%%
if ~exist('data_type','var')
  data_type='bid';
end
%{
while (~(strcmp(data_type(1:3),'tra')...
        |strcmp(data_type(1:3),'nca')...
        |strcmp(data_type(1:3),'mar')...
        |strcmp(data_type(1:3),'lej')...
        |strcmp(data_type(1:3),'hyg')...
        |strcmp(data_type(1:3),'imm')...
        |strcmp(data_type(1:3),'sim')...
        |strcmp(data_type(1:3),'gui'))),
%}
switch data_type(1:3)
    case {'tra','nca','mar','lej','hyg','imm','sim','gui'}, %OK
    otherwise %not accepted data_type
  display('You have not specified the data type. You must set :')
  display('            data_type="ncar"     for EISCAT data in NCAR format')
  display('            data_type="guisdap"  for EISCAT data in GUISDAP format')
  display('            data_type="lejeune"  for EISCAT data in Lejeune format')
  display('            data_type="mars"     for data coming from mars model')
  display('            data_type="transcar" for data coming from transcar model')
  display('            data_type="hygraph"  for EISCAT data in Lejeune format')
  display('')
  %data_type=input('Enter the value for data_type : ','s');
  %data_type=data_type(1:3);
  error('set data_type and try again')
% pour transcar il faut definir file
% pour guisdap  il faut definir path_expr

end %switch
%%
while ~exist('data_path','var')
   if strcmp(data_type(1:3),'gui') || strcmp(data_type(1:3),'nca')
      % Debut Modif DA/03 Aout 1998
      if ~exist('outpath','var')
        data_path=input('Enter the Data Path  [ Hit "Return" if local ] : ','s');
      else
	data_path=outpath;
      end
      % Fin Modif DA/03 Aout 1998
   else
      data_path=input('Enter the Data Path : ','s');
   end
end
%%

%pc=computer;
if ~isempty(data_path)
  len=length(data_path);char=data_path(len);
  data_path(len+1) = '/';
%   if char~=':' & char~='/' & char~='\'
%     switch pc(1:2)
%      case 'PC'
%        data_path(len+1)='\';
%      case 'MA'
%        data_path(len+1)=':';
%      otherwise
%        data_path(len+1)='/';
%     end
%   end
end


clear('len','char')

switch(data_type(1:3))
    case 'nca'
        chemin_donnees=dir([data_path '*.*dtst']); 
        if isempty(chemin_donnees) || chemin_donnees(1).bytes==0
            display(['Sorry: path ',data_path,...
                  ' not exist or is empty of *.ndtst files'])
            clear('data_path')
        end %if
    case 'gui'
        chemin_donnees=dir([data_path 'filelist.dat']);
        if isempty(chemin_donnees) || chemin_donnees(1).bytes==0
            display(['Sorry: path ',data_path,...
                  ' not exist or is empty of filelist.dat files'])
            clear('data_path')
        end
end
%% read data from TRANSCAR

switch data_type(1:3)
    case 'tra', read_tra;
    otherwise, error(['case ',data_type,' not yet handled.'])
end %switch
%%      

end_hour=time(nb_records);
str_hour=time(1);

%if nb_records>1
%  init_figures

%  default_limits

  req_pos=1;
%%
  set_position
%%
%  display_figures
%  if nb_pos > 4
%     if ~remote
%        plotscan
%     end
%  end
%end
