%OPEN THE FILE
fid=fopen(isfile,'r');



%HOW MANY ALTITUDE POINTS?
sline=fgets(fid); sline=fgets(fid);   %ff past labels
isnz=0;
while length(sline)>=6 && ~strcmp(sline(1:6),'Sondre')
    sline=fgets(fid);
    isnz=isnz+1;
end
isnz=isnz-1              %there is a space between records
fseek(fid,0,'bof');      %rewind the file



%LOOP OVER INDIVIDUAL DATA RECORDS READING AND ORGANIZING THEM
isit=1;
while ~feof(fid)
    %status update
    if mod(isit,50) == 0
        fprintf('\nISREAD.M --> reading record #:  %d',isit);
    end
    
    %ff past labels
    sline=fgets(fid); sline=fgets(fid);

    %loop over individual lines of data and check for goodness
    for isiz=1:isnz
        sline=fgets(fid);
        numline=str2num(sline);

        %there is some missing or bad data
        if isempty(numline)
            minds=strfind(sline,'missing');
            badinds=strfind(sline,'knownBad');

            for im=1:length(minds)
                sline(minds(im):minds(im)+6)='0000000';
            end
            for ib=1:length(badinds)
                sline(badinds(ib):badinds(ib)+7)='00000000';
            end
            numline=str2num(sline);
        end

        %organize the data in this line
        istime(isit)=numline(1);
        isz(isiz)=numline(4);
        
        isne(isiz,isit)=numline(6);
        isdne(isiz,isit)=numline(7);
        iste(isiz,isit)=numline(8);
        isdte(isiz,isit)=numline(9);
        isti(isiz,isit)=numline(10);
        isdti(isiz,isit)=numline(11);
        isvi(isiz,isit)=numline(12);
        isdvi(isiz,isit)=numline(13);
    end
    sline=fgets(fid);       %trailing blank line
    
    %time index counter
    isit=isit+1;
end
fprintf('\n');



%CLOSE FILE
fclose(fid);



%SAVE DATA
save(['./',isfile,'.mat'],'is*');