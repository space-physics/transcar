function read_kinmod(datadir,kinfile,outfile)

%
%function read_kinmod(datadir,kinfile,outfile)
%

%OPEN FILE
fid=fopen([datadir,'/',kinfile],'r');



%ARRAY INITIALIZATION
nenref=0; tset=[];                                  %for keeping up with kinetic code energy grid adjustments vs. time
day=[]; sec=[]; dipangle=[]; nalt=[]; nen=[];       %header info (data sizes, etc.)
zop=[]; no1d=[]; no1s=[]; noii2p=[]; nn2a3=[]; po3p3p=[]; po3p5p=[]; p1ng=[]; pmein=[]; p2pg=[]; p1pg=[];        %population densities (nx) or rates (px)
e=[]; fluxdown=[];
%clear e* fluxdown*;



%MAIN LOOP FOR READING IN TIME SERIES
k=1;
while(~feof(fid))
    if(mod(k,50)==0)
        fprintf('KINREAD.M --> Reading record #:  %d\n',k);
    end

    head=fscanf(fid,'%e',5);                                                                                               %read in header for this time    
    if(feof(fid))
        break;
    end
    day=[day,head(1)]; sec=[sec,head(2)]; dipangle=[dipangle,90-head(3)]; nalt=[nalt,head(4)]; nen=[nen,head(5)];             %store header info in appropriate vars.
    
    data=fscanf(fid,'%e',[11 nalt(k)])';           %read in main data portion

    %organize into appropriately name arrays
    zop=[zop,data(:,1)]; 
    no1d=[no1d,data(:,2)]; 
    no1s=[no1s,data(:,3)]; 
    noii2p=[noii2p,data(:,4)]; 
    nn2a3=[nn2a3,data(:,5)]; 
    po3p3p=[po3p3p,data(:,6)]; 
    po3p5p=[po3p5p,data(:,7)];
    p1ng=[p1ng,data(:,8)]; 
    pmein=[pmein,data(:,9)];
    p2pg=[p2pg,data(:,10)]; 
    p1pg=[p1pg,data(:,11)];

    %read in precipitation information
    if(nen(k) ~= nenref)
        tset=[tset,k];
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

eval(['e',int2str(tsetl),'=e;'])
eval(['fluxdown',int2str(tsetl),'=fluxdown;']);
clear e; clear fluxdown;
fclose(fid);



%FIX TIME VARIABLE SO THAT IT DOESN'T WRAP AROUND AT MIDNIGHT
timeop=sec/3600;
ith=max(find(timeop>=23));
iws=find(timeop(ith:length(timeop))<23);
timeop(ith+iws-1)=timeop(ith+iws-1)+24;

end
