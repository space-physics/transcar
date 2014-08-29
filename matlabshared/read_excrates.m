function read_excrates(datadir,kinfile,outfile)
%
% Michael Hirsch: You may wonder why the while/fscanf paradigm is
% used--it's because the original data from emissions.dat is written in
% 1960's punch-card 80 column style--data wraps around from line to line in
% a way not usable by textscan!

%% OPEN FILE
kinFN = [datadir,'/dir.output/',kinfile];
display(['reading:  ',kinFN])
fid=fopen(kinFN,'r');
%% ARRAY INITIALIZATION
%for keeping up with kinetic code energy grid adjustments vs. time
nenref=0; tset=[];
%header info (data sizes, etc.)
%day=[]; sec=[]; dipangle=[]; nalt=[]; nen=[];
%population densities (nx) or rates (px)
%zop=[]; no1d=[]; no1s=[]; noii2p=[]; nn2a3=[]; po3p3p=[]; po3p5p=[]; p1ng=[]; pmein=[]; p2pg=[]; p1pg=[];
e=[]; fluxdown=[];
%clear e* fluxdown*;

%% MAIN LOOP FOR READING IN TIME SERIES
k=1;
while(~feof(fid))
%     if(mod(k,200)==0)
%         fprintf('READ_EXCRATES.M --> Reading record #:  %d\n',k);
%     end

    head=fscanf(fid,'%e',5);                                                                                               %read in header for this time    
    if(feof(fid))
        break;
    end

    %store header info in appropriate vars.
    day(k)     = head(1); 
    sec(k)     = head(2); 
    dipangle(k)= 90-head(3); 
    nalt(k)    = head(4); 
    nen(k)     = head(5);  

    %read in main data portion
    data=fscanf(fid,'%e',[11 nalt(k)])';    

    %organize into appropriately named arrays
    zop(:,k)   = data(:,1); 
    no1d(:,k)  = data(:,2); 
    no1s(:,k)  = data(:,3); 
    noii2p(:,k)= data(:,4); 
    nn2a3(:,k) = data(:,5); 
    po3p3p(:,k)= data(:,6); 
    po3p5p(:,k)= data(:,7);
    p1ng(:,k)  = data(:,8); 
    pmein(:,k) = data(:,9);
    p2pg(:,k)  = data(:,10); 
    p1pg(:,k)  = data(:,11);

    %read in precipitation information
    if(nen(k) ~= nenref)
        tset = [tset,k];
        tsetl=length(tset);
        nenref=nen(k);
        if(tsetl>1)
            eval(['e',int2str(tsetl-1),'=e;'])
            eval(['fluxdown',int2str(tsetl-1),'=fluxdown;']);
        end
        e=[]; fluxdown=[];
    end
    data=fscanf(fid,'%e',[2 nen(k)])';
    e=[e,data(:,1)];
    fluxdown=[fluxdown,data(:,2)];

    k=k+1;
end

eval(['e',int2str(tsetl),'=e;']) %FIXME: yuck
eval(['fluxdown',int2str(tsetl),'=fluxdown;']); %FIXME: yuck
clear('e','fluxdown')
fclose(fid);



%FIX TIME VARIABLE SO THAT IT DOESN'T WRAP AROUND AT MIDNIGHT
timeop=sec/3600;
ith=find(timeop>=23, 1, 'last' );
iws=find(timeop(ith:length(timeop))<23);
timeop(ith+iws-1)=timeop(ith+iws-1)+24; %#ok<NASGU>



%SAVE THE EXCITATION RATES TO AN OUTPUT .MAT FILE
outFN = [datadir,'/',outfile];
display(['READ_EXCRATES.M --> Saving calculations to: ',outFN])
save(outFN,'n*','p*','timeop','zop','dipangle','fluxdown*','e*')

end %function
