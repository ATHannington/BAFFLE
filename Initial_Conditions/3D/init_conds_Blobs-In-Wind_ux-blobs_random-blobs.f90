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
		Real*8 :: delx,dely,delz					!Ranges in x, y and z for the maximum range limit
		Real*8 :: d							!Thickness of zero velocity layer

		real*8,external::isotherm_p

		Real*8 :: exp_blob,Dmax,Rhomax,RhoMin
		Real*8 :: x0,y0,z0,rhoblob,dblob

		Real*8 :: xarg,yarg,zarg,dist
		Integer :: i
		Real*8 :: p_BackGround,rho_BackGround				!Background Pressure & density for Adiabatic and Two Phase Isothermal Setups

!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)	


	delx=xrangeMax-xrangeMin
	dely=yrangeMax-yrangeMin
	delz=zrangeMax-zrangeMin


!	Set thickness of zero velocity layer:
	d=0.375d0

!	Set background phase Pressure and Density:
	p_BackGround=rhominset*(CsIso**2)
	rho_BackGround=rhominset

!	Set e blob to 2 to achieve circular regions:
	exp_blob=2
!	Set max diameter of blob:
	Dmax=DmaxSet
!	Set max density of blob:
	rhomax=RhomaxSet
!	Set Background Density : 
	RhoMin=RhominSet

	If ((RandBlobBool.eqv..true.).AND.(BlobSwitchBool.eqv..TRUE.)) then
!	Allocate arrays the dimension of number of blobs being used:
		allocate (x0_blob_array(NumBlobs))
		allocate (y0_blob_array(NumBlobs))	
		allocate (z0_blob_array(NumBlobs))	
		allocate (d_blob_array(NumBlobs))	
		allocate (rho_blob_array(NumBlobs))	

		call rand_blobs(d,dmax,rhomax)
	endif
!	Switch off any further calls to generate random blob regions:
	RandBlobBool=.false.

!-------------------------------------------------------------------------------
!	Setup Background Atmosphere
!-------------------------------------------------------------------------------

	rho = RhoMin

!	Set Isothermal Pressure
	If (IsoBool.eqv..TRUE.) then
		p = isotherm_p(rho)
	else
		p=p_BackGround
	endif


	If ((x>=0.d0).and.(x<=0.d0+d)) then
		If (IsoBool.eqv..TRUE.) then
			ux=+1.5d0*DSQRT((p_BackGround/rho_BackGround))
		else if (IsoBool.eqv..FALSE.) then
			ux=+1.5d0*DSQRT(gam*(p/rho))
		endif
		uy=0.d0
		uz=0.d0
	else if ((x>=delx-d).and.(x<=delx)) then
		If (IsoBool.eqv..TRUE.) then
			ux=-1.5d0*DSQRT((p_BackGround/rho_BackGround))
		else if (IsoBool.eqv..FALSE.) then
			ux=-1.5d0*DSQRT(gam*(p/rho))
		endif
		uy=0.d0
		uz=0.d0
	else
		ux=0.d0
		uy=0.d0
		uz=0.d0
	endif

!-------------------------------------------------------------------------------
!	Setup Blobs based on randomly generated conditions:
!-------------------------------------------------------------------------------
	IF (BlobSwitchBool.eqv..TRUE.) then
		Do while (i<=NumBlobs)
	!	Get data from blob arrays previously randomly generated:
			x0 = x0_blob_array(i)
			y0 = y0_blob_array(i)
			z0 = z0_blob_array(i) 
			dblob = d_blob_array(i)
			rhoblob = rho_blob_array(i)

	!	Dist(x-vec,o-vec) = [((x_x-o_x)/(0.5d0*d_x))**e+((x_y-o_y)/(0.5d0*d_y))**e+((x_x-o_z)/(0.5d0*d_z))**e]**(1/e)
	!	If Dist<=1 point is inside region
			zarg=((z-z0)/(0.5d0*dblob))**exp_blob
			yarg = ((y-y0)/(0.5d0*dblob))**exp_blob
			xarg = ((x-x0)/(0.5d0*dblob))**exp_blob
			dist = (zarg+yarg+xarg)**(1.d0/exp_blob)

			If ((x0.gt.(dblob+0.02d0)).and.(x0.lt.(delx-(dblob+0.02d0))).and.(dist<=1.d0)) then
				rho = rhoblob
				uy= 0.d0
				uz=0.d0
				If (IsoBool.eqv..TRUE.) then
					p = isotherm_p(rho)
				else
					p=p_BackGround
				endif
			endif
			i=i+1
		enddo
	endif

end subroutine init_conds


!===============================================================================


!===============================================================================
!#####		RANDOMLY GENERATE DENSE COLD BLOB REGIONS		########
!===============================================================================

subroutine rand_blobs(d,dmax,rhomax)
	use constants								!Tells program to use global constants
	implicit none

		Real*8,intent(in) :: d,dmax,rhomax
		Real*8 :: dblob,rho,exp_blob,x0,y0,z0
		Real*8 :: rand_num
		Real*8 :: delx,dely,delz					!Ranges in x, y and z for the maximum range limit
		Integer :: i


	delx=xrangeMax-xrangeMin
	dely=yrangeMax-yrangeMin
	delz=zrangeMax-zrangeMin

	call srand(RandSeed)

	i=1

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!	Generate Blobs:
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	Do while (i<=NumBlobs)
		101 CONTINUE
		Print'(A15,1x,I3.3,1x,A2,1x,I3.3)',"Generating Blob",i,"of",NumBlobs
		
!	Generate a random number:
		rand_num = rand(0)
!-------------------------------------------------------------------------------
!	Setup Diameter:
!-------------------------------------------------------------------------------
		If (rand_num.gt.0.d0) then
!	Set diameter of blob:
			dblob = rand_num*dmax
		else 
!	The loop will restart to call a new random number if the diameter is zero
			GOTO 101
		endif
!-------------------------------------------------------------------------------
!	Set Centre of region:
!-------------------------------------------------------------------------------
		102 CONTINUE
		rand_num = rand(0)
!		If the x cooridinate lies outside of the wind regions a new x coord will be calculated:
!		Also checks that the x coord is not too close to the boundary:
		If ( (((rand_num*delx).gt.(dblob+0.02d0)).and.((rand_num*delx).lt.(d)))& 
		.or.(((rand_num*delx).lt.(delx-(dblob+0.02d0))).and.((rand_num*delx).gt.(delx-d)))) then
				x0 = rand_num*delx
		else
			GOTO 102
		endif


		103 CONTINUE
		rand_num = rand(0)
!		If the y cooridinate lies too close to the bounaries:
		If (((rand_num*dely).gt.(dblob+0.02d0)).and.((rand_num*dely).lt.(dely-(dblob+0.02d0)))) then
			y0 = rand_num*dely
		else
			GOTO 103
		endif


		104 CONTINUE
		rand_num = rand(0)
!		If the z cooridinate lies too close to the bounaries:
		If (((rand_num*delz).gt.(dblob+0.02d0)).and.((rand_num*delz).lt.(delz-(dblob+0.02d0)))) then
			z0 = rand_num*delz
		else
			GOTO 104
		endif

!-------------------------------------------------------------------------------
!	Set Density of Region:
!-------------------------------------------------------------------------------
		If (TwoPhaseBool .eqv. .FALSE.) then
			105 CONTINUE
			rand_num = rand(0)
	!		If the density is non-zero and above critical density for two phase, set density as fraction of maximum density:
			If ((rand_num.gt.0.d0).and.((rand_num*RhoMax).ge.RhoCrit)) then
				rho = rand_num*rhomax
			else
				GOTO 105
			endif
		else if (TwoPhaseBool .eqv. .TRUE.) then
			rho=rhomax
		endif
!-------------------------------------------------------------------------------
!	Set Variables into arrays:
!-------------------------------------------------------------------------------
!	Allocate above variables to blob arrays:
	x0_blob_array(i) = x0
	y0_blob_array(i) = y0
	z0_blob_array(i) = z0
	d_blob_array(i) = dblob
	rho_blob_array(i) = rho

		i=i+1
	enddo 

end subroutine rand_blobs
