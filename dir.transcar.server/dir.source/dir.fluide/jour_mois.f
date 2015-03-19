       subroutine jour_mois(ian,ijour,imois)
Cf2py intent(in) ian
Cf2py intent(out) ijour
Cf2py intent(out) imois
Cf2py intent(out) iday

       integer m1(12),m2(12)

       data m1/31,59,90,120,151,181,212,243,273,304,334,365/
       data m2/31,60,91,121,152,182,213,244,274,305,335,366/

       ijour=mod(ian,1000)
       nan=(ian-ijour)/1000
       njour=mod(nan,4)

	do 10 i=1,12

	  if (njour.eq.0) then
	    if (ijour.gt.m2(i)) goto 10
            if (i.eq.1) then
            ijour=ijour
            else   
	    ijour=ijour-m2(i-1)
            endif
            imois=i  
	    goto 20
	  else
	    if (ijour.gt.m1(i)) goto 10
            if (i.eq.1) then
            ijour=ijour
            else   
	    ijour=ijour-m1(i-1)
            endif
            imois=i 
	    goto 20
	  endif

10	continue
20	continue

	return

	entry mois_jour(ian,iday)

	nan=ian/10000
	njour=mod(nan,4)
	iday=ian-nan*10000
	im=iday/100
	iday=iday-100*im
	if (im.gt.1) iday=iday+m1(im-1)
	if (njour.eq.0.and.im.ge.3) iday=iday+1

	return
	end
