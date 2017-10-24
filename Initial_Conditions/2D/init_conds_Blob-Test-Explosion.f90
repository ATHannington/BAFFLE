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
		Real*8 :: r_blob, r_explsn					!Radius in x-y plane of blob, and explosion point
		Real*8 :: delx,dely,delz					!Ranges in x, y and z for the maximum range limit
		Real*8 :: lowbound,midbound,upbound				!Bounds for ranges in initial conditions
		Real*8 :: x0_blob,y0_blob,x0_explsn,y0_explsn			!x and y coordinates of centre of high pressure region of blob, and explosion point
		Real*8 :: dmax_blob,dmax_explsn					!Maximum radius of high pressure region of blob, and explosion

		Real*8 :: dist_blob,dist_explsn
		Real*8 :: xarg, yarg
		Integer :: exp_blob,exp_explsn

!	Cells format: 	cells(i,j,k,1:23) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z,grad_ux,grad_uy,grad_uz, &
!						& grad_p_x,grad_p_y,grad_p_z	/)	

!	Set x0_blob, y0_blob to be on diagonal of system:
	x0_blob=0.57d0
	y0_blob=0.57d0

!	Set max diameter of blob:
	dmax_blob=0.10d0

!	Set e blob to 2 to achieve circular region:
	exp_blob=2

!	Dist(x-vec,o-vec) = [((x_x-o_x)/(0.5d0*d_x))**e+((x_y-o_y)/(0.5d0*d_y))**e+((x_x-o_z)/(0.5d0*d_z))**z]**(1/e)
!	If Dist<=1 point is inside region
	xarg=((x-x0_blob)/(0.5d0*dmax_blob))**exp_blob
	yarg = ((y-y0_blob)/(0.5d0*dmax_blob))**exp_blob
	dist_blob = (xarg+yarg)**(1.d0/exp_blob)
	

!	Set radius of blob:
!	r_blob = Sqrt((x-x0_blob)**2+(y-y0_blob)**2)!+z**2)
!	r_blob = Sqrt(x**2+y**2)	



!	Set x0_explsn, y0_explsn to be at centre of system:
	x0_explsn=0.50d0
	y0_explsn=0.50d0

!	Set max radius of high pressure region of explosion:
	dmax_explsn=0.09591663046d0

!	Set e blob to 2 to achieve circular region:
	exp_explsn=2

!	Dist(x-vec,o-vec) = [((x_x-o_x)/(0.5d0*d_x))**e+((x_y-o_y)/(0.5d0*d_y))**e+((x_x-o_z)/(0.5d0*d_z))**z]**(1/e)
!	If Dist<=1 point is inside region
	xarg=((x-x0_explsn)/(0.5d0*dmax_explsn))**exp_explsn
	yarg = ((y-y0_explsn)/(0.5d0*dmax_explsn))**exp_explsn
	dist_explsn = (xarg+yarg)**(1.d0/exp_explsn)


!	Set radius of explosion:
!	r_explsn = Sqrt((x-x0_explsn)**2+(y-y0_explsn)**2)!+z**2)
!	r_blob = Sqrt(x**2+y**2)


!	Set axis lengths and calculate lower bound, mid-pont, and upper bound:
	delx=xrangeMax-xrangeMin
	dely=yrangeMax-yrangeMin
	delz=zrangeMax-zrangeMin

	lowbound=0.d0
	midbound=0.5d0
	upbound=Sqrt(delx**2+dely**2+delz**2)

!Sets values of pressure, x-velocity, and density according to starting conditions for 2D Sod Shock
	If (dist_blob<=1.d0) then
		rho=1.d0
		uy=0.d0
		uz=0.d0
		p=0.1d0
		ux=0.d0
	else if (dist_explsn<=1.d0) then
		rho=0.5d0
		ux=0.d0
		uy=0.d0
		uz=0.d0
		p=2.d0
	else
		rho=0.5d0
		uy=0.d0
		uz=0.d0
		p=0.1d0
		ux=0.d0
	endif

end subroutine init_conds
!===============================================================================
