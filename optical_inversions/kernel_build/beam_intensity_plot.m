%GENERATE A HISTOGRAM OF BEAM ENERGIES, MAINLY USED TO CHECK THAT THE
%SIMULATIONS RAN WITH THE CORRECT INPUT.  YOU MUST RUN PROCESS_BEAMS ON THIS 
%SET OF BEAM SIMULATIONS BEFORE RUNNING THIS FILE.
basedir='~/simulations/TRANSCAR/beamsims/';
%datelab='01Mar2011v2';
datelab='10Dec2007';
load([basedir,datelab,'/ids.mat']);


for l=1:length(energies)-1
    if(energies(l)==floor(energies(l)))
      beamenergy=[num2str(energies(l)),'.'];
    else
      beamenergy=num2str(energies(l));
    end
    fprintf('LOADING DATA FOR BEAM:  %f eV \n',energies(l));
    load([basedir,datelab,'/beam',beamenergy,'/excrates.mat'],'e*','fluxdown*','time*');

    its=find(timeop>onsettime);
    E(:,l)=e1(:,its(1));
    phi(:,l)=fluxdown1(:,its(1));
end

figure;
loglog(E,phi,'LineWidth',2)
set(gca,'FontSize',16);
xlabel('energy (eV)');
ylabel('downward hemisphere diff. num. flux (cm^{-2} s^{-1} eV^{-1})')
axis([50,20e3,1e3,1e10]);
