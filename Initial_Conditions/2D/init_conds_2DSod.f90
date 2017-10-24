!===============================================================================
!
!*** 				INITIAL CONDITIONS			     ***
!
!-------------------------------------------------------------------------------
!Subroutine for initial conditions functions for density, velocity, and pressure
!& as a function of x.
!-------------------------------------------------------------------------------
subroutine init_conds(x,y,z,rho,ux,uy,uz,p)			
	use constants								!Tells program to use global constants
	implicit none

		Real*8, intent(in) :: x,y,z
		Real*8, intent(out) :: rho,ux,uy,uz,p
		Real*8 :: r_xy							!Radius in x-y plane
		Real*8 :: delx,dely,delz					!Ranges in x, y and z for the maximum range limit
		Real*8 :: lowbound,midbound,upbound				!Bounds for ranges in initial conditions
		Real*8 :: x0,y0							!x and y coordinates of centre of high pressure region
		Real*8 :: rmax							!Maximum radius of high pressure region

!	Cells format: 	cells(i,j,k,1:23) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z,grad_ux,grad_uy,grad_uz, &
!						& grad_p_x,grad_p_y,grad_p_z	/)	

!	uy=0.d0
	uz=0.d0

!	Set x0, y0 to be centre of system:
	x0=0.5d0
	y0=0.5d0

!	Set max radius of high pressure region:
	rmax=0.50d0/2.d0

!	Set radius:
	r_xy = Sqrt((x-x0)**2+(y-y0)**2)!+z**2)
!	r_xy = Sqrt(x**2+y**2)	

!	Set axis lengths and calculate lower bound, mid-pont, and upper bound:
	delx=xrangeMax-xrangeMin
	dely=yrangeMax-yrangeMin
	delz=zrangeMax-zrangeMin

	lowbound=0.d0
	midbound=0.5d0
	upbound=Sqrt(delx**2+dely**2+delz**2)

!Sets values of pressure, x-velocity, and density according to starting conditions for 2D Sod Shock
	If ((r_xy>=0.d0).and.(r_xy<rmax)) then
		rho=1.d0
		ux=0.d0
		uy=0.d0
		uz=0.d0
		p=1.d0
	else if ((r_xy>=rmax)) then
		rho=0.125d0
		ux=0.d0
		uy=0.d0
		uz=0.d0
		p=0.1d0
	else
		Print*,"Invalid coordinate! Program will terminate!"		!Exits program if x lies outside of 1 to 0 range
		STOP	
	endif

end subroutine init_conds
!===============================================================================
