%function [alt,en,t,fhemd,fhemu]=read_diffnumflux(filename)
%
%This function reads in a TRANSCAR output file and organizes the simulation
%results into upward and downard flux distributions.
%
%   INPUT:
%       filename - name of TRANSCAR output file containing differential
%           number flux vs. altitude, time, and energy
%
%   OUTPUTS:
%       alt - altitude grid for differential number flux (1D array, km)
%       en - energy grid (1D array, eV)
%       t - Universal time (1D array, s)
%       fhemd - avg. downward pitch angle differential number flux (3D array, cm^-2 s^-1 eV^-1 sR^-1)
%       fhemu - avg. upward pitch angle differential number flux (3D array, cm^-2 s^-1 eV^-1 sR^-1)
%
%The dimension ordering of output differential number flux arrays is
%altitude (row dim.), energy (column dim.), time (3rd dim.).


function [alt,en,t,fhemd,fhemu]=read_diffnumflux(filename)

fid=fopen(filename,'r');


t=[];
fhemd=[];
fhemu=[];
l=1;
while ~feof(fid)
    %SOME USER OUTPUT
    if mod(l,10)==0
        fprintf('READ_DIFFNUMFLUX.M --> Reading record #:  %d\n',l);
    end
    
    %DATA SIZE INFO
    metadata=fscanf(fid,'%f %d %d',3);
    if isempty(metadata)
       break; 
    end
    t=[t,metadata(1)];
    nalt=metadata(2);
    nen=metadata(3);


    
    %DISTRIBUTION ENERGY GRID
    en=fscanf(fid,'%f',[1 nen]);
    fluxup=zeros(nalt,nen);
    fluxdown=zeros(nalt,nen);
    
    
    
    %LOOP OVER ALTITUDES AND READ IN UP AND DOWN DIST. FN. AT EACH GRIDPOINT
    for k=1:nalt
        %ALTITUDE GRID
        alt(k)=fscanf(fid,'%f',1);

        %DOWNWARD AND UPWARD DIFFERENTIAL NUMBER FLUX (CM-2 S-1 EV-1 SR-1)
        %VS ENERGY AT THIS ALTITUDE
        fluxes=fscanf(fid,'%e',[2 nen]);
        fluxdown(k,:)=fluxes(1,:);
        fluxup(k,:)=fluxes(2,:);
    end
    
    %REORDER ALTITUDE DIMENSION TO BE FROM LOW TO HIGH
    alt=flipud(alt(:));
    fluxdown=flipud(fluxdown);
    fluxup=flipud(fluxup);
    
    
    
    %PUT INTO AN OUTPUT DATA CUBE
    fhemd=cat(3,fhemd,fluxdown);
    fhemu=cat(3,fhemu,fluxup);
    l=l+1;
end

fclose(fid);

end