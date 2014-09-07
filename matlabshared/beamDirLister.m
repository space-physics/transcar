function idlist = beamDirLister(basedir,energies)
%Michael Hirsch 2014

nEnergies = length(energies);

for iEnergy=1:nEnergies
    idlist{iEnergy,1}=[basedir,'/beam',num2str(energies(iEnergy)) ];
end %for

%FIXME: do we really need to make floating points end in decimal? Does the
%Fortran code require this?
addDecInd = cellfun(@isempty,regexp(idlist,'.*/beam\d*\.','match'));

[decPts{1:nEnergies,1}] = deal(''); %initialize
[decPts{addDecInd,1}]   = deal('.'); %stuff decimals
idlist = strcat(idlist,decPts); %actually add the trailing decimal point where needed

end %function