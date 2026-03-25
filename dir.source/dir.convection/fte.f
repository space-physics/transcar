	subroutine fte(tu,latmag_ref,tmag_ref,latmag,tmag,Bmag,VM,WM,
     &   		dpot,dEx,dEy)



	data pi/3.14/

	real x,y
	real gamma,Dx_gamma,Dy_gamma
	real f1,f2,Dx_f1,Dy_f1,Dx_f2,Dy_f2
	real Pthese
	real R,Dx_r2,Dy_r2
	real latmag,tmag,latmag_ref,tmag_ref
 	real tu,m,n,B0,V0,beta0

	logical in

	real tmag2dist,latmag2dist

	real Ex_in,Ex_out,dEx,Ey_in,Ey_out,dEy
	real pot_in,pot_out,dpot

	real W_in,W_out,W,v_in,v_out,v
	real Vx_in,Vx_out,Vy_in,Vy_out,Vx,Vy

c	ellipse's half axis toward magnetic est

	m = 350.

c	ellipse's half axis toward magnetic north

	n = 100.

c tmag2dist converts 1hMLT into dist. in km

	tmag2dist   = 2.*pi*6370.*cos(latmag)/24.

c latmag2dist is the dist.in km corresponding to 1degre of lat.

	latmag2dist = 2.*pi*6370./360.

	x = (tmag-tmag_ref)*tmag2dist
	y = (latmag-latmagref)*latmag2dist

	in = ((x/n)**2+(y/m)**2).le.1

	B0 = Bmag
	V0 = sqrt(VM**2+WM**2)


	Pthese =  (x**2-y**2+m**2-n**2)

	gamma  =  0.5*atan2(2*x*y,Pthese)
	if (x.lt.0) gamma=gamma+pi
       	Dx_gamma = 0.5*(2*y*Pthese-4*(x**2)*y)/Pthese**2/
     &	 	(1+(2*x*y/Pthese)**2)
	Dy_gamma = 0.5*(2*x*Pthese+4*x*y**2)/Pthese**2/
     &		(1+(2*x*y/Pthese)**2)

	R  =  (Pthese**2+4*x**2*y**2)
	Dx_r2 = 0.25*(4*x*Pthese+8*x*y**2)*R**(-0.75)
	Dy_r2 = 0.25*(-4*y*Pthese+8*x**2*y)*R**(-0.75)

	f1 = R**(0.25)*cos(gamma)-x
	f2 = R**(0.25)*sin(gamma)-y
	Dx_f1 = -sin(gamma)*Dx_gamma*R**(0.25)+cos(gamma)*Dx_r2-1
	Dx_f2 = cos(gamma)*Dx_gamma*R**(0.25)+sin(gamma)*Dx_r2
	Dy_f1 = -sin(gamma)*Dy_gamma*R**(0.25)+cos(gamma)*Dy_r2
	Dy_f2 = cos(gamma)*Dy_gamma*R**(0.25)+sin(gamma)*Dy_r2-1


c Plasma velocity

	Vx_out = -V0/(m-n)*(n*sin(beta0)*Dy_f1+m*cos(beta0)*Dy_f2)
	Vy_out = V0/(m-n)*(n*sin(beta0)*Dx_f1+m*cos(beta0)*Dx_f2)

	Vx_in = V0*cos(beta0)
	Vy_in = V0*sin(beta0)



c Velocity potential


	W_in = -V0*(cos(beta0)*x+sin(beta0)*y)
	W_out = -V0/(m-n)*(n*f2*sin(beta0)-m*f1*cos(beta0))


c Electric potential

	dpot = -B0*W

c Streamline function

	v_in = V0*(cos(beta0)*y-sin(beta0)*x)
	v_out = -V0/(m-n)*(n*f1*sin(beta0)+m*f2*cos(beta0))


c Electric field (x- and y-components)


	Ex_in = -V0*B0*sin(beta0)
	Ey_in = V0*B0*cos(beta0)

	Ex_out = -V0*B0/(m-n)*(n*Dx_f1*sin(beta0)+m*Dx_f2*cos(beta0))
	Ex_out = -V0*B0/(m-n)*(n*Dy_f1*sin(beta0)+m*Dy_f2*cos(beta0))



	if (in) then
		Vx = Vx_in
		Vy = Vy_in

		W = W_in
		v = v_in

		dpot = pot_in

		dEx = dEx_in
		dEy = dEy_in

     	else
		Vx = Vx_out
		Vy = Vy_out

		W = W_out
		v = v_out

		dpot = pot_out

		dEx = dEx_out
		dEy = dEy_out

	end if

	return
	end
