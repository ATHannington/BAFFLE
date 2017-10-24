!-------------------------------------------------------------------------------
module constants								!Module declaring all global constants	

!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)

!	Define constants:
	Integer, save :: NcXr,NcX,NcYr,NcY,NcZr,NcZ,Ncellsg=2,Ncells		!Number of Cells: real, in x,y and z,(NcXr etc.) and their total with ghost counterparts (NcX etc.);
										!ghost (for implementation of boundary conditions); total
	Real*8, allocatable,dimension(:,:,:,:) :: cells				!Cells of simulation array over which the variables are sampled
	Real*8 :: xrangeMax,xrangeMin						!Min and max of x-axis range of test volume
	Real*8 :: yrangeMax,yrangeMin						!Min and max of y-axis range of test volume
	Real*8 :: zrangeMax,zrangeMin						!Min and max of z-axis range of test volume
	Real*8 :: gam								!Gamma, the ratio of specific heats Cp/Cv
	Real*8 :: delta_t,delt_default=1.d-4						!Time step for Time Integration, and default to prevent zeros being used.
	Real*8 :: t								!Time up to which the integrator will proceed
	Real*8 :: CellVol							!Cell Volume
	Real*8 :: C_CFL								!C_CFL: Courant-Friedrichs-Lewy number such that 0<C_CFL<1
	Real*8 :: Cs								!Sound speed
	Real*8 :: xArea, yArea, zArea						!Surface area for x,y,z direction faces of cuboid.
	Integer :: mn,qxn,qyn,qzn,En,rhon,uxn,uyn,uzn,pn,Vn,xn,yn,zn, &		!Integers to hold the position of each in the d-dimension (cells(i,j,k,d)) of cells array.
 		 & gdxn,gdyn,gdzn, &
		 & guxxn,guyxn,guzxn, &
		 & guxyn,guyyn,guzyn, &
		 & guxzn,guyzn,guzzn, &
	 	 & gpxn,gpyn,gpzn

	Real*8 :: Pi=DAcos(-1.d0)

!	Set save file output directory:
	character(len=100):: SavePath = "/home/alpha/ath4/summer2017/HydroCode/D_TEST/"		
!***PLEASE NOTE***: A path string of greater than 100 characters will NOT be tolerated within the program unless you edit the hard-coded '(A100,~)' format identifiers in the write_data subroutine in hydrocode-main, to '(Ax,~)' where x>100
!	Set File format of outputs
	character(len=100) :: SaveFileFormat = ".dat"				


!	Switch questions on or off. Reverts to following defaults if off.
!	Off is useful for running the program as a background job.
	Logical :: Qbool=.FALSE.						!Toggle to switch runtime questions on or off:	true => questions off	false => questions on

!	Settings for switching between accuracy orders:
	Logical :: orderBool							!Integer to hold what order of evaluation has been selected. FALSE => First Order TRUE=> Second Order
	Logical :: orderBooldefault= .true.					!Second Order: FALSE => First Order TRUE=> Second Order

!	Some system constants defaults:
	Real*8 ::NcXrdefault=100.d0,NcYrdefault=1.d0,NcZrdefault=1.d0
	Real*8 :: tdefault=0.2d0
	Real*8 :: gamDefault=5.0d0/3.d0
	Real*8 :: CFLDefault=0.4d0

!	Frame writing settings:
	Logical :: writeFramesbool						!Integer to hold whether to write frames or just write at end. true => write frames ; false => write at end only
	Logical :: wFramesDefault= .TRUE.					! true => Write Frames	 false => Don't write Frames
	Logical :: delt_from_frames_bool = .TRUE.				!TRUE => Adjust timesteps with number of frames (can dramatically slow program, but will enable more visualisation at 											small final times). This also ensures the number of frames selected 
										!FALSE => Do not adjust the timestep from number of frames irrespective of whether frames are selected.
!	Set Max number of Frames to be written:
	Integer :: MaxFrames = 100


!	Isothermal Settings:
	Logical :: IsoBool							!" True=> Isothermal Setup ON	False => Adiabatic Setup ON"
	Logical :: IsoBoolDefault= .FALSE.					!" True=> Isothermal Setup ON	False => Adiabatic Setup ON"
	Real*8 :: CsIso								!Isothermal Global Sound-Speed
	Real*8 :: CsIsoDefault= DSQRT(0.1d0/0.5d0)				!Default Isothermal Global sounds speed. Edit this based on CsIso = SQRT(p_desired / RhoMaxSet)
	LOGICAL :: TwoPhaseBool = .FALSE.					!Toggles Two Phase Isothermal Setup. TRUE => Two Phase On	FALSE => Two Phase Off
	Real*8 :: IsoDeltaTScaler = 1.d0/DSQRT(5.d0/3.d0)			!Scales Isothermal timestep to be equivalent to that of Adiabtic case.		

!	Random Number Generator Seed:
	Integer*4 :: RandSeed
	Integer*4 :: RandSeedDefault=1234 

!	Set maximum and minimum density, and diameter of blobs:
!	Diameter:
	Real*8 :: DmaxSet=0.2d0
!	Density:
	Real*8 :: RhomaxSet=2.0d0
	Real*8 :: RhominSet=0.5d0						
	Real*8 :: RhoCrit							!Value Halfway between Rho min and max that sets the minimum value of Rho that Blobs can take

!	Call to Random Blobs Bool:
	Logical :: RandBlobBool=.TRUE.						!This switch allows the random blob regions to be defined once, on the first call to init_conds, and no other time after 											that. 	***MUST ALWAYS BE TRUE***
	LOGICAL :: BlobSwitchBool = .TRUE.					!This toggles switches whether blobs are to be generated or not.	TRUE => Blobs Generated		FALSE => BLobs 											NOT Generated

	Integer :: NumBlobs = 25						!Number of blobs to be generated in the random blob generator
	Real*8, allocatable,dimension (:) :: x0_blob_array,y0_blob_array, &	!Arrays to hold the parameters of the blobs globally such that they can be used on every call of init_conds
				& z0_blob_array, d_blob_array, rho_blob_array




end module constants

!-------------------------------------------------------------------------------
