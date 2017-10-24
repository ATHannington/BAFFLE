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
	Integer :: orderBool							!Integer to hold what order of evaluation has been selected.
	Real*8 :: Pi=DAcos(-1.d0)
	integer :: writeFramesbool						!Integer to hold whether to write frames or just write at end. =1 => write frames ; =0 => write at end only


	integer :: Qbool=1							!Toggle to switch runtime questions on or off:	0 => questions off	1 => questions on
!	Deafult values if questions are turned off:
	Integer :: orderBooldefault=1 						!Second Order **ON** :: 0 => First Order 1=> Second Order
	Real*8 ::NcXrdefault=100.d0,NcYrdefault=100.d0,NcZrdefault=1.d0
	Real*8 :: tdefault=0.2d0
	Real*8 :: gamDefault=5.0d0/3.d0
	Real*8 :: CFLDefault=0.4d0
	Integer :: wFramesDefault=1						! 1 => Write Frames	 0 => Don't write Frames

!	Isothermal Settings:
	Integer :: IsoBool							!" 1=> Isothermal Setup ON	0 => Adiabatic Setup ON"
	Integer :: IsoBoolDefault=1						!" 1=> Isothermal Setup ON	0 => Adiabatic Setup ON"
	Real*8 :: CsIso								!Isothermal Global Sound-Speed
	Real*8 :: CsIsoDefault=1.d0						!Default Isothermal Global sounds speed
end module constants

!-------------------------------------------------------------------------------
