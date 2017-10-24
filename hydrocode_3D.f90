!Andrew Hannington								
!ath4
!Created: 29 May 2017
!Last edited: 22 June 2017
!3D Hydro-code, for solving of simple cases such as Sod Shock Tube

!-------------------------------------------------------------------------------
!*******************************************************************************
!-------------------------------------------------------------------------------

program hydrocode								!Main program

	use constants								!Tells program to use global constants
	implicit none								!No implicit definitions used

!	Define variables:

!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)		

!	Initialise constants:
	call init_constants							!Calls subroutine to define global constant values

!	Initialise constants in Riemann Solver:
	call initialise_gammas

!	Initialise cells array:
	call cell_setup

!	Apply boundary conditions:
	call bounds

!	Performs time integral using primitive variable conversion and riemann
!	solver, from time t=0 to time t=0.2
	call time_int

!	Write Data at end only if user has selected not to print FRAMES:
	If (writeFramesBool.eqv..false.) then
		call write_data(0,t)
	endif

end program hydrocode
!-------------------------------------------------------------------------------
!*******************************************************************************
!-------------------------------------------------------------------------------

!Subroutine of initial constants:
!-------------------------------------------------------------------------------
subroutine init_constants							

	use constants
	implicit none
		character(len=1) :: timeQ,gamQ,NcQ,charTest,orderQ,cflQ,isoQ, &
				  & CsQ,writeQ
		Integer :: error,charLen				!Takes error code for read
		Character(len=4):: error_loc				!Error Location character string· Identifies error location for error subroutine

		Real*8 :: t_Suggested						!Suggested End-Time
		Real*8 :: x_Dist						!Distance proportion for Winds to travel in x

!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)			

!Integers to hold the position of each in the d-dimension (cells(i,j,k,d)) of cells array:
 		mn=1
		qxn=2
		qyn=3
		qzn=4
		En=5
		rhon=6
		uxn=7
		uyn=8
		uzn=9
		pn=10
		Vn=11
		xn=12
		yn=13
		zn=14
		gdxn=15
		gdyn=16
		gdzn=17
		guxxn=18
		guyxn=19
		guzxn=20
		guxyn=21
		guyyn=22
		guzyn=23
		guxzn=24
		guyzn=25
		guzzn=26
		gpxn=27
		gpyn=28
		gpzn=29

!	Sets up the error location identifier for this subroutine:
	error_loc="cons"

	If(qBool.eqv..true.) then
	!	The following section allows the user to select second or first order evaluation:
		Print*,"Select first (f) or second (s) order evaluation:"
		Read*, orderQ

	!		Set global boolean to contain this choice:
		If ((orderQ.eq."f").or.(orderQ.eq."F")) then
			orderBool=.false.
		else if ((orderQ.eq."s").or.(orderQ.eq."S")) then
			orderBool=.true.
		else
			Print*,"Order Select Boolean Not Selected! Default Set:",orderBooldefault
			Print*,"FALSE => First Order	TRUE => Second Order"
			orderBool=orderBoolDefault
		endif

	else if (qBool.eqv..false.) then
		orderBool=orderBooldefault
	endif



	If(qBool.eqv..true.) then
	!	Following section allows user to adjust number of spatial cells the range is divided into

		Print*,"Would you like to adjust number of spatial cells? (y/n)"
		Read*, NcQ
	
		If ((NcQ=="y").or.(NcQ=="Y")) then
			Print*,"Enter number of cells in x-direction:"
			Read*,NcXr

			Print*,"Enter number of cells in y-direction:"
			Read*,NcYr

			Print*,"Enter number of cells in z-direction:"
			Read*,NcZr

	!error_sub(IO_error,error_loc,NcXrdum,NcYrdum,NcZrdum,gamdum,m,E,V,p1,p2,d1,d2,i,j,k,time)
	!	Here we are only checking the dimensions of real cells selected by user, so all other values
	!	set to dummy, non-error producing values.		
			call error_sub(0,error_loc,NcXr,NcYr,NcZr,0.5d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1,1,1,1.d0)

		else
			Print*,"Default number of cells selected. Number of Cells in each direction:"
			Print*,"x=",NcXrdefault,"y=",NcYrdefault,"z=",NcZrdefault
			NcXr=NcXrdefault
			NcYr=NcYrdefault
			NcZr=NcZrdefault
		endif

	else if (qBool.eqv..false.) then
		NcXr=NcXrdefault
		NcYr=NcYrdefault
		NcZr=NcZrdefault
	endif


	NcX=NcXr+Ncellsg
	NcY=NcYr+Ncellsg
	NcZ=NcZr+Ncellsg

	Ncells=NcX*NcY*NcZ							!Number of cells is equal to number of real cells plus number of ghost boundary cells
	allocate (cells(NcX,NcY,NcZ,29))					!Allocate the number of cells the space is divided into (Ncells) to the number of 
										!& rows to be used (cells), and the number of variables to be stored (23) to the number of columns.

!	Sets ranges of x, y and z coordinates				
	xrangeMax=1.d0
	xrangeMin=0.d0

	yrangeMax=1.d0
	yrangeMin=0.d0

	zrangeMax=1.d0
	zrangeMin=0.d0


	If(qBool.eqv..true.) then
	!	User enters gamma values. A default of 5/3 for strong shocks used if none selected. 
	!&	Program aborts if invalid entries or unphysical values selected.
		Print*,
		Print*,"Do you want to enter specific heats constant, gamma? (y/n)"
		Read*,gamQ

		If ((gamQ=="y").or.(gamQ=="Y")) then
			Print*,"Enter Gamma value (5/3>gamma>1):" 
			Read(*,*,IOSTAT=error) gam
	
	!error_sub(IO_error,error_loc,NcXrdum,NcYrdum,NcZrdum,gamdum,m,E,V,p1,p2,d1,d2,i,j,k,time)
	!		Here we are checking the error code for the validity of the user input
	!		(non-numerical values will be flagged up in the error subroutine),
	!		and the value of gamma that has been entered. If this is non-physical
	!		this will also be flagged by the error subroutine.	
			call error_sub(error,error_loc,NcXr,NcYr,NcZr,gam,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1,1,1,t)

		else
			Print*,"Default Gamma Set:",gamDefault
			gam=gamDefault
		endif

	else if(qBool.eqv..false.) then
		gam=gamDefault
	endif



!	Suggests end-time for supersonic wind and blobs calculation:

	If(qBool.eqv..true.) then
	!	User inputs end time for integration. A default is selected otherwise and system
	!&	aborts if invalid, or unphysical entries are given.
		Print*,"Do you want to adjust time integration upper limit?(y/n)"
		Read*, timeQ
		
		If ((timeQ=="y").or.(timeQ=="Y")) then
			Print*,"Okay! Enter time (t) in seconds:"		
			read(*,*,IOSTAT=error) t					!IOSTAT returns 0 for normal operation, and non-zero for errors. Allows for program to be aborted with invalid entry.

	!error_sub(IO_error,error_loc,NcXrdum,NcYrdum,NcZrdum,gamdum,m,E,V,p1,p2,d1,d2,i,j,k,time)
	!		Here we are checking the error code for the validity of the user input
	!		(non-numerical values will be flagged up in the error subroutine),
	!		and the value of time that has been entered. If this is non-physical
	!		this will also be flagged by the error subroutine.	
			call error_sub(error,error_loc,NcXr,NcYr,NcZr,0.5d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1,1,1,t)

		else
			Print*,"Final Time set to default:",tdefault
	!		DefaultEnd time for integrator
			t=tdefault !0.2d0
		endif

	else if(qBool.eqv..false.) then

		t=tdefault

	endif

	If(qBool.eqv..true.) then

		Print*,"Do you want to adjust CFL number for integration accuracy control?(y/n)"
		Read*, cflQ
		If ((cflQ=="y").or.(cflQ=="Y")) then
			Print*,"Okay! CFL number (0<=CFL<=0.5):"		
			read(*,*,IOSTAT=error) C_CFL
			If (error.ne.0) then
				Print*,"Invalid Entry! Program aborted!"
				STOP
			else if (C_CFL.lt.0.d0) then
				Print*,"Invalid Entry (C_CFL<0)! Program aborted!"
				STOP
			else if (C_CFL.gt.0.5d0) then
				Print*,"Warning! Unstable CFL Selected! Program may fail!"
			endif
		else if ((cflQ=="n").or.(cflQ=="N")) then
			Print*,"Default CFL number set! C_CFL=",CFLDefault
		!	Set C_CFL: Courant-Friedrichs-Lewy number such that 0<C_CFL<1
			C_CFL=CFLDefault
		else
			Print*,"Default C_CFL set:",CFLDefault
			C_CFL=CFLDefault
		endif

	else if(qBool.eqv..false.) then
	
		C_CFL=CFLDefault

	endif


!	Question to select isothermal or adiabatic settings
	If(qBool.eqv..true.) then

		Print*,"Do you want Isothermal (i) or Adiabatic (a) settings?"
		Read*, isoQ
		If ((isoQ=="i").or.(isoQ=="I")) then
			IsoBool = .true.
			Print*,"ISOTHERMAL Selected!"

			Print*,
			Print*,"Do you want to adjust Global Sound-Speed (Cs)? (y/n)"
			Print*,"Current Default Cs =",CsIsoDefault

			Read*, CsQ
			If ((CsQ=="Y").or.(CsQ=="Y")) then
				Print*,"Enter Global Sound-Speed! Cs="
				Read*,CsIso
			else if ((CsQ=="n").or.(CsQ=="N")) then
				Print*,"Default Global Sounds-Speed (Cs) Selected!"
				CsIso = CsIsoDefault

				Print*,"Cs=",CsIso
			else
				Print*,"Default Cs set:",CsIsoDefault
				CsIso = CsIsoDefault
			endif


		else if ((isoQ=="a").or.(isoQ=="A")) then
			IsoBool=.false.
			Print*,"ADIABATIC Selected!"
		else
			Print*,"Default set:",IsoBoolDefault
			Print*," TRUE=> Isothermal Setup ON	FALSE => Adiabatic Setup ON"
			IsoBool=IsoBoolDefault
			
			If(IsoBool.eqv..true.) then
				Print*,"ISOTHERMAL: Default Global Cs set:",CsIsoDefault
				CsIso = CsIsoDefault
			endif
		endif

	else if(qBool.eqv..false.) then
	
		IsoBool=IsoBoolDefault
		CsIso = CsIsoDefault

	endif

!	Set Initial timestep:
	delta_t=delt_default

! 	Set dummy sound speed to prevent machine precision being used in error:
	Cs=1.d0


!	Cell volume:
	CellVol = ((xrangeMax-xrangeMin)/NcXr)*((yrangeMax-yrangeMin)/NcYr)*((zrangeMax-zrangeMin)/NcZr)

!	Set surface areas of faces
	xArea=((zrangeMax-zrangeMin)/NcZr)*((yrangeMax-yrangeMin)/NcYr)
	yArea=((xrangeMax-xrangeMin)/NcYr)*((zrangeMax-zrangeMin)/NcZr)
	zArea=((xrangeMax-xrangeMin)/NcXr)*((yrangeMax-yrangeMin)/NcYr)

!	Set Blobs Random Number seed:
	RandSeed=RandSeedDefault

!	Set Writing multiple outputs option:
	If(qBool.eqv..true.) then
		Print*,"Do you want to write multiple outputs of snapshots in time (FRAMES)? (y/n)"
		Read*, writeQ

		If ((writeQ=="y").or.(writeQ=="Y")) then
			writeFramesBool= .true.
		else if ((writeQ=="n").or.(writeQ=="N")) then
			writeFramesBool= .false.
		else
			Print*,"Invalid Entry! Write Frames Default Setting used:",wFramesDefault
			Print*,"False => No Frames Written	True => Frames Will Be Written"
			writeFramesBool=wFramesDefault
		endif
	else if(qBool.eqv..false.) then
		writeFramesBool=wFramesDefault
	endif

!	Print Write Frames Setting Selected:
	Print*,"Write Frames setting selected:",writeFramesBool
	Print*,"TRUE => Write Frames	FALSE => Do NOT Write Frames"
	Print*,

	Print*,
	Print*,"Maximum Number of Data Outputs (FRAMES):", MaxFrames

	RhoCrit = 1.0d0 * ((RhoMaxSet-RhoMinSet)/2.d0)

end subroutine init_constants

!-------------------------------------------------------------------------------
!Subroutine of setting up cell starting values
!-------------------------------------------------------------------------------
subroutine cell_setup
	use constants								!Tells program to use global constants
	implicit none

		Real*8 :: m,qx,qy,qz,E,rho,ux,uy,uz,p
!		m:mass, q:fluid momentum components, E:Total Energy, 
!&		rho:Density, u:fluid velocity components, p:pressure
		Real*8 :: Midx, Midy, Midz, V					!Midpoint of cell in cartesian coordinates and cell volume
		Real*8 :: delx,dely,delz					!System dimensions
		Integer :: i,j,k						

!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)		
	V= CellVol

!	Temporarily adjust system dimensions
	delx=(xrangeMax-xrangeMin)*(1.d0+(2.d0/(NcXr)))
	dely=(yrangeMax-yrangeMin)*(1.d0+(2.d0/(NcYr)))
	delz=(zrangeMax-zrangeMin)*(1.d0+(2.d0/(NcZr)))

!	Spatial counter:
	i=1									
	j=1
	k=1

!	Assign starting values from init_conds, and mid-points to cells array:
	Do while (k<=NcZ)
		j=1
		Do while (j<=NcY)
			i=1
			Do while (i<=NcX)

				If (i==1) then
					Midx=(delx/(NcX))*(-0.5d0)
					Midy=(dely/(NcY))*((j-1)-0.5d0)
					Midz=(delz/(NcZ))*((k-1)-0.5d0)!Sets the middle of the cell to be at the middle of the kth fraction of the z
				else if (j==1) then
					Midy=(dely/(NcY))*(-0.5d0)
					Midx=(delx/(NcX))*((i-1)-0.5d0)!Sets the middle of the cell to be at the middle of the ith fraction of the x range
					Midz=(delz/(NcZ))*((k-1)-0.5d0)!Sets the middle of the cell to be at the middle of the kth fraction of the z
				else if (k==1) then
					Midz=(delz/(NcZ))*(-0.5d0)
					Midx=(delx/(NcX))*((i-1)-0.5d0)!Sets the middle of the cell to be at the middle of the ith fraction of the x range
					Midy=(dely/(NcY))*((j-1)-0.5d0)!Sets the middle of the cell to be at the middle of the jth fraction of the y range

				else
					Midx=(delx/(NcX))*((i-1)-0.5d0)!Sets the middle of the cell to be at the middle of the ith fraction of the x range
					Midy=(dely/(NcY))*((j-1)-0.5d0)!Sets the middle of the cell to be at the middle of the jth fraction of the y range
					Midz=(delz/(NcZ))*((k-1)-0.5d0)!Sets the middle of the cell to be at the middle of the kth fraction of the z range			
				endif

!			Sets values of zero for boundary cells, as otherwise they will be left with undefined values due to init_conds not applying
!			outside of the real cell range:			
				If ((i==1).or.(i==NcX).or.(j==1).or.(j==NcY).or.(k==1).or.(k==NcZ)) then
					cells(i,j,k,1:29) = (/0.d0,qx,qy,qz,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,V,Midx,Midy,Midz, &
							&	0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
							&	0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
				else
		!			Call initial conditions to evaluate density, pressure and velocity
		!&			at midpoint of each cell:
					Call init_conds(Midx,Midy,Midz,rho,ux,uy,uz,p)

		!			Initial conditions for mass, and energy:
					m = V*rho
					E = V*( (p/(gam-1.d0)) + 0.5d0*rho*((ux**2)+(uy**2)+(uz**2)) )

		!Set initial momenta:
					qx= ux*m
					qy= uy*m
					qz= uz*m

		!			Update cells:
					cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
							&	0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
							&	0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)	

				endif

				i=i+1
			enddo
			j=j+1
		enddo
		k=k+1
	enddo

end subroutine cell_setup

!===============================================================================
!#####		Compact subroutine for time integration			########
!===============================================================================
subroutine time_int
	use constants								!Tells program to use global constants
	implicit none

		real*8 :: time							!Counter to track the time progression. Is passed to subroutines for error message purposes.	
		Real*8 :: twopercent,timedum_prcnt				!value of two percent of final time, and dummy to track value up to two percent
		Integer :: frames						!Number of frames to be written	
		Real*8 :: print_time						!time at which write is performed
		Real*8 :: print_time_tracker					!Tracker of next time at which write is to be performed
		Real*8 :: delt_tol						!********Tolerance of MINIMUM VALUE OF DELTA_T**********
		Real*8 :: percent						!Percentage complete tracker
		Integer :: i,j,k
		Integer :: frame_num						!Frame (output) number 
		Real*8 :: Etot							!Total Energy variable, an output from variable conversion for error checking purposes
		LOGICAL :: vcBool						!Variable Conversion Bool. Switched off by t=0 write frame to skip double variable conversion in that case.
		LOGICAL :: deltBool						!Delta_t Boolean. Carries whether the delta_t of timesteps is adjusted by the number of frames
		Integer :: tsteps_taken
!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)	

!	Total number of frames to be written:
	Frames = maxFrames

!	Set time at which data is written:
	print_time = t/dble(Frames)
	print_time_tracker = print_time

!********Tolerance of MINIMUM VALUE OF DELTA_T**********
	delt_tol = 1.d-5

!	Updated time 
	time = 0.d0
	timedum_prcnt=0.d0		

!	Calculate two percent of final time:
	twopercent = t*0.02d0


!	If user has selected .TRUE. in Constants, then the timestep will be adjusted by the number
!	of frames written iff frames are selected to be written. Else, the program does not adjust
!	the timestep irrespective of whether the frames are to be written.
	If (delt_from_frames_bool.eqv..TRUE.) then
		deltBool = writeFramesBool
	else
		deltBool = .FALSE.
	endif

!	Initialise boolean for variable conversion
	vcBool= .TRUE.

!	Initialise frame numbers:	
	frame_num = 0

	If (writeFramesBool.eqv..true.) then
	!	Write initial conditions
		call var_conv(time,Etot)
!Print*,"t,del_t,Etot:",time,delta_t,Etot
		call write_data(frame_num,time)
!	Set variable conversion to false such that variable conversion is not performed twice in first run
		vcBool=.FALSE.
	endif

!	Set up tracker for counting how many timesteps have been taken.
!	Set to 1 because the multiplier for the timestep needs to be non-zero.
!	This tracker is brought back in line with number of steps taken by logic
!	following time update.
	tsteps_taken = 1



!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!#####			TIME INTEGRATION				########
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
	Do while (time<t)
		
		If (vcBool.eqv..TRUE.) then
!			Convert conserved to primitive variables:
			call var_conv(time,Etot)
		endif

!-------------------------------------------------------------------------------
! Delta_t Adjusting Logic
!-------------------------------------------------------------------------------
!
!	To ensure only the final time is reached, and no further, this condition sets the final timestep equal
!	to the difference between current and end times:
		If (delta_t.gt.(t-time)) then
			delta_t = t-time
		endif
		
!	If the option to adjust delta_t by the max frame number has been selected
!	(internally checked against whether the write frames has been selected in
!	the first place [see logic above]), the timestep is adjusted backwards to
!	being only a maximum of the write time intervals if too large, and is
!	adjusted smaller (as timesteps can NEVER be increased) to match to the
!	write time when the next timestep will take it beyond the write time.
!					--
!	There is a maxval condition that sets the lowest value delta_t can take.
!	This is to ensure against infinite loops when time==Print_time_tracker
!	i.e. current simulation time is equal to next write time.
		If (deltBool .eqv. .true.) then		
			If (delta_t.gt.print_time) then
				delta_t = print_time
			else if (delta_t.lt.print_time) then
				If ((time+delta_t).gt.Print_time_tracker) then
					delta_t = Print_time_tracker - time	
				endif

				delta_t = maxval((/delta_t,delt_tol/))
			endif
		endif
!-------------------------------------------------------------------------------

!Print*,"t,del_t,Etot:",time,delta_t,Etot

!		Pass cells through Riemann Solver:
		call cells_to_riemann(time)

!		Update time counter:
		time = time + delta_t
		
!		Update number of timesteps taken, ignoring the first step to bring tracker back
!		 in track with time.
		If (time>0.d0) then
			tsteps_taken = tsteps_taken + 1
		endif

!		Calculate percentage complete
		percent = (time/t)*100.d0

!		Update time dum tracker
		timedum_prcnt = timedum_prcnt + delta_t


		If(timedum_prcnt>=twopercent) then

			Print*,
			Print'(A21,I3,A1)',"Percentage complete= ",int(ceiling(percent)),"%"
		!Reset percentage tracking counter to print further 2%'s
			timedum_prcnt=0.d0

			Print'(A2,1x,ES10.2E3)',"t=",time
		endif


		If (writeFramesBool.eqv..true.) then
			If (time.ge.print_time_tracker) then
!				Update Frame Number:
				frame_num = frame_num + 1

	!			Write the data for given frame:
				call write_data(frame_num,time)		
				print_time_tracker = print_time_tracker + print_time
			endif
		endif


!	Set variable conversion bool back to true such that all following timesteps perform variable conversion:
		vcBool=.TRUE.

	enddo
!-------------------------------------------------------------------------------

	Print*,
	Print*,"***TOTAL*** Number of Data Outputs (FRAMES):",(frame_num+1)
	Print*,"(Max Frames + 1 t=0 Frame)"

end subroutine time_int

!-------------------------------------------------------------------------------
!Subroutine to write data to output files
!-------------------------------------------------------------------------------
subroutine write_data(frame_num,time)
	use constants								!Tells program to use global constants
	implicit none

		Integer, intent(in) :: frame_num
		Real*8,intent(in):: time
		integer :: i,j,k
		character(len=:),allocatable :: left_final
		character(len=:),allocatable :: filename
		character(len=:),allocatable :: left,mid,right
		character(len=200) :: tmp1,tmp2,tmp3,tmp4
		Integer :: val
		LOGICAL :: LengthAdjustBool

!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)	


!	write(1,*) (/"/__x_/","/__y_/","/__z_/","/__m_/","/_qx_/","/_qy_/","/_qz_/",&
!			"/__E_/","/rho_/","/_ux_/","/_uy_/","/_uz_/", &
!			"/__p_/","/__V_/"/)
!

!	Setup left and right of filename:
	Left="Output_"
	Mid = "_t="
	Right= SaveFileFormat
!-------------------------------------------------------------------------------
!	Generate SavePath & File name string:
!-------------------------------------------------------------------------------
	left_final = trim(adjustl(SavePath)) // trim(adjustl(left))

	write(tmp1,'(A100,I5.5)') trim(adjustl(left_final)),frame_num
	
!	tmp2 = trim(adjustl(tmp1)) // trim(adjustl(Mid))

!	write(tmp3,'(A100,ES9.2E3)') trim(adjustl(tmp2)),time

	tmp4 = trim(adjustl(tmp1)) // trim(adjustl(Right))

	filename =  trim(adjustl(tmp4))
!-------------------------------------------------------------------------------
	
	Open(1,file=trim(adjustl(filename)))

!	Insert time to each file:
	write(1,'(A2,1x,ES15.6E3)') "t=",time

!	Insert format header to each data file:
	write(1,*) (/"   x  ","   y  ","   z  ","   m  ","  qx  ","  qy  ","  qz  ",&
			"   E  "," rho  ","  ux  ","  uy  ","  uz  ", &
			"   p  ","   V  "/)

!	i=2 so data written skips ghost cell. Similarly for ncells-1 in do loop
	i=2
	j=2
	k=2

	Print'(A13,1x,I10,1x,A2,1x,I10)',"Writing Frame",frame_num,"of",MaxFrames

	Do while(k<=NcX-1)
		j=2
		Do while (j<=NcY-1)!NcY-1)
			i=2
			Do while(i<=NcX-1)	
				write(1,'(20(ES28.18E4,5x))') (/cells(i,j,k,xn),cells(i,j,k,yn),cells(i,j,k,zn),Cells(i,j,k,1:11)/)
				i=i+1
			enddo
			j=j+1
		enddo
		k=k+1
	enddo

	close(1)

	Print*,"Finished Writing Frame!"

end subroutine write_data
!-------------------------------------------------------------------------------
! Variable Conversion routine. Derives the primitive variables of density, velocity
! and pressure from the conserved variables of mass, volume, momentum, energy and
! fluid velocity.
!-------------------------------------------------------------------------------
subroutine var_conv(time,Etot)
	use constants								!Tells program to use global constants
	implicit none
		Real*8, intent(in) :: time
		Real*8:: ux,uy,uz						!Velocities are edited and passed out for time step stability criterion
		Real*8 :: rho,p							!pressure and density to be brought in, used in calculation, and passed out to time int, then to timestep, in order to 
										!calculate soundspeed
		Character(len=4) :: error_loc					!Error Location character string· Identifies error location for error subroutine	
		Real*8 :: m,qx,qy,qz,E
!		m:mass, q:fluid momentum components, E:Total Energy, 
!&		rho:Density, u:fluid velocity components, p:pressure
		Real*8 :: Midx,Midy,Midz, V					!Midpoint of cell and cell volume
		Integer :: i,j,k
		Real*8 :: grad_d_x,grad_d_y,grad_d_z, &
			  grad_ux_x,grad_uy_x,grad_uz_x, &
			  grad_ux_y,grad_uy_y,grad_uz_y, &
			  grad_ux_z,grad_uy_z,grad_uz_z, &
			  grad_p_x,grad_p_y,grad_p_z

		Real*8 ,intent(out):: Etot					!Total Energy of the system
		Real*8 :: Etotdum
		real*8,external::isotherm_p					!Isothermal Pressure

!	Declare location indicator for error subroutine:
	error_loc="varC"

!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)	
	Etot=0.d0


	i=1
	j=1
	k=1

	Do while (k<=NcZ)
		j=1
		Do while(j<=NcY)
			i=1
			Do while (i<=NcX) 
			
		!		Variable conversion:

				m = cells(i,j,k,mn)
				qx = cells(i,j,k,qxn)
				qy = cells(i,j,k,qyn)
				qz = cells(i,j,k,qzn)
				E = cells(i,j,k,En)
				rho = cells(i,j,k,rhon)
				ux = cells(i,j,k,uxn)
				uy = cells(i,j,k,uyn)
				uz = cells(i,j,k,uzn)
				p = cells(i,j,k,pn)
				V = cells(i,j,k,Vn)
				Midx = cells(i,j,k,xn)
				Midy = cells(i,j,k,yn)
				Midz = cells(i,j,k,zn)
				grad_d_x= cells(i,j,k,gdxn)
				grad_d_y=cells(i,j,k,gdyn)
				grad_d_z=cells(i,j,k,gdzn)
				grad_ux_x=cells(i,j,k,guxxn)
				grad_uy_x=cells(i,j,k,guyxn)
				grad_uz_x=cells(i,j,k,guzxn)
				grad_ux_y=cells(i,j,k,guxyn)
				grad_uy_y=cells(i,j,k,guyyn)
				grad_uz_y=cells(i,j,k,guzyn)
				grad_ux_z=cells(i,j,k,guxzn)
				grad_uy_z=cells(i,j,k,guyzn)
				grad_uz_z=cells(i,j,k,guzzn)
				grad_p_x=cells(i,j,k,gpxn)
				grad_p_y=cells(i,j,k,gpyn)
				grad_p_z=cells(i,j,k,gpzn)

		!		rho = m/V:
				rho=m/v

		!		u=q/m:
				ux = qx/m
				uy = qy/m
				uz = qz/m

		!		p = (gamma - 1)*(E/V - 1/2 * rho * u^2)
				p = (gam-1.d0)*( (E/V)-( 0.5d0*rho*((ux**2)+(uy**2)+(uz**2)) ) )


!			Set Isothermal Pressure:
				If (IsoBool.eqv..TRUE.) then
					p = isotherm_p(rho)

				endif

				Etotdum = Etotdum + E

		!		Checks for zero or negative pressures to prevent NaNs:

		!Checks for zero or negative pressures, masses, volumes and energies to prevent NaNs:
!		error_sub(IO_error,error_loc,NcXrdum,NcYrdum,NcZrdum,gamdum,m,E,V,p1,p2,d1,d2,i,j,k,time)
				call error_sub(0,error_loc,NcXr,NcYr,NcZr,gam,m,E,V,p,1.d0,1.d0,1.d0,i,j,k,time)

				cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
						& grad_d_x,grad_d_y,grad_d_z, &
						& grad_ux_x,grad_uy_x,grad_uz_x, &
						& grad_ux_y,grad_uy_y,grad_uz_y, &
						& grad_ux_z,grad_uy_z,grad_uz_z, &
						& grad_p_x,grad_p_y,grad_p_z /)	

!				Calculate stable timestep:
				call timestep(i,j,k,ux, uy, uz,p,rho,time)

				i=i+1
			enddo
			j=j+1
		enddo
		k=k+1		
	enddo

!	Apply boundary conditions:
	call bounds

!	If second order has been selected, gradients will be calculated.
	If (orderBool.eqv..TRUE.) then
	!	Calculate New gradients:
		Call get_gradients(time)

	!	Apply boundary conditions:
		call bounds
	endif


	Etot=Etotdum

end subroutine var_conv
!-------------------------------------------------------------------------------
! Time Step stability subroutine. Ensures cells can only speak to nearest
!& neighbours following: delta_t<C_CFL *(delta_x/v_signal)
! Delta_t: timestep
! C_CFL: Courant-Friedrichs-Lewy number such that 0<C_CFL<1
! delta_x: size of single cell
! v_signal: maximum signal of velocity at that time V_signal=V+Cs
! In N dimensions: delta_t < C_CFL * Sigma[i=1->N](delta_x_i/v_i_signal)
!
! This subroutine selects the smallest of all stable timesteps from the total
! system, and sets that as the sytem-wide timestep, using indices from Var_conv.
!-------------------------------------------------------------------------------
subroutine timestep(i,j,k,ux, uy, uz,p,rho,time)
	use constants								!Tells program to use global constants
	implicit none
		Integer, intent(in) :: i,j,k					!Cell counters coming in from variable conversion
		Real*8, intent(in) :: ux,uy,uz					!velocities are passed in from time integration conversion to calculate max timestep
		Real*8, intent(in) :: p,rho					!pressure and density passed in to calculate sound speed
		Real*8, intent(in) :: time	
		Real*8 ::delx,dely,delz,delmin					!Size of cell in x, y, and z directions
		Real*8 :: deltmax						!Maximum value of delta_t
		Real*8 :: reduced_delta_t
		Real*8 :: umag							!Magnitude of velocity

!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)	

!	Calculate dimensions of each cell:
	delx=((xrangeMax-xrangeMin)/(NcXr))
	dely=((yrangeMax-yrangeMin)/(NcYr))
	delz=((zrangeMax-zrangeMin)/(NcZr))

	delmin=minval((/delx,dely,delz/))

!	Calculate sound speed for given cell:
	If (IsoBool .eqv. .FALSE.) then
		cs = DSQRT(gam*(p/rho))	
	else
		Cs = DSQRT(p/rho)
	endif

	umag=DSQRT(ux**2+uy**2+uz**2)

!	Calculate maximum time step for given cell
	deltmax = C_CFL*(delmin/(umag+Cs))

	If(IsoBool.eqv..FALSE.)then
		reduced_delta_t = deltmax
	else if (IsoBool.eqv..TRUE.) then
		reduced_delta_t = IsoDeltaTScaler * deltmax
	endif

!	For first cell only, set new minimum timestep:
	If ((i==1).and.(j==1).and.(k==1)) then
		delta_t=reduced_delta_t
!	If newly calucated timestep is less than previous smallest timestep, set timestep
!	to new, smaller value.
	else if (deltmax.lt.delta_t) then
		delta_t=reduced_delta_t
!	If new timestep is no smaller than previous smallest, timestep is not adjusted:
	else

	endif

end subroutine timestep

!-------------------------------------------------------------------------------
!Subroutine to pass cells through Riemann solver (credit to K. Wood, B. Vandenbroucke,
! E. Toro). 
!-------------------------------------------------------------------------------
!* @brief Solve the Riemann problem with the given left and right state.
!*
!* @param d1 Left state density.
!* @param u1 Left state fluid velocity.
!* @param p1 Left state pressure.
!* @param d2 Right state density.
!* @param u2 Right state fluid velocity.
!* @param p2 Right state pressure.
!* @param gg Adiabatic index of the gas.
!* @param ds Variable to store the resulting density in.
!* @param us VReal*8, intent(in): timeariable to store the resulting fluid velocity in.
!* @param ps Variable to store the resulting pressure in.
!* @param uflag Flag variable used to signal if the solution was sampled from the
!* left state (-1) or the right state (1). A vacuum state (0) is currently not
!* supported by this version of the Riemann solver.
!* @param x_over_t Point (in velocity space) where we want to sample the
!* solution. In a finite volume scheme you usually want to set this value to 0.
!------------------------------------------------------------------------------

subroutine cells_to_riemann(time)
	use constants								!Tells program to use global constants
	implicit none
		Real*8, intent(in) :: time
		Character(len=4) :: error_loc				!Error Location character string· Identifies error location for error subroutine
		Real*8, dimension(3) :: Avec					!Surface Area vector
		Real*8 :: magA							!magnitude of surface area vector, A.
		Real*8, dimension(3) :: norm					!Surface Normal vector
		Character(len=1) :: xyz						!Character passed in from cells to Riemann that tells this subroutine whether we are looking at x, y, or z direction.
		Integer :: i,j,k						!Integers from Cells to that identify which cells the calculation needs data for
		Integer :: uflag						!Indicates whether the left or right cell other velocity components should be used for solution velocity. E.g. if 											!on x-axis analysis us from Riemann is Uvec,x-component, and if uflag=-1 u_y=u_cell,y, u_z=u_cell,z,
										!else if uflag=+1 u_y=u_cell+1,y, u_z=u_cell+1,z
		Real*8 :: x_over_t						!" @param x_over_t Point (in velocity space) where we want to sample the
										! solution. In a finite volume scheme you usually want to set this value to 0. "
										! As such we will set this value to be zero here.
		Real*8 :: u1x,u1y,u1z,d1,p1,u2x,u2y,u2z,d2,p2,gg,ds,ps
		Real*8 :: us,usx,usy,usz
		Real*8 :: Ax,Ay,Az						!Surface normal. A=1 if left position lower than right positiom
										!A=-1 if left position higher in x/y/z than right position
		Real*8, dimension(3):: u1vec,u2vec,usvec			!Vectors to contain the components of velocity for Left hand velocity, RH vel, and solution vel.

		Real*8 :: dL,pL,dR,pR,uL,uR
		Real*8, dimension(3) :: uLvec,uRvec

!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)	
	gg= gam
	x_over_t=0.d0
	Avec=(/0.d0,0.d0,0.d0/)
	
!	Feed cells through Riemann solver, with 1->Ncells being right to left.
!	Ncells-1 to stop iteration going outside number of cells in use. 


!	Define left and right pressures, densities, and velocities for the Riemann solver, imposing reflective boundary conditions.
	


!	 cell_select (xyz,uflag,i,j,k,d1,u1,p1,d2,u2,p2)

	i=1
	j=1
	k=1

	Do while (k<=NcZ-1)
		j=1
		Do while (j<=NcY-1)
			i=1
			Do while (i<=NcX-1)
!------------------------------------------------------------------------------
!			Select x direction:				
				xyz='x'

!		Setup Unit Area Vector A-Vec. Only x-direction non-zero for x-axis flux:
				Ay=0.d0
				Az=0.d0
		!		Check left to right x positions:
				If (cells(i,j,k,xn).lt.cells(i+1,j,k,xn)) then
					Ax=1.d0*xArea
				else if (cells(i,j,k,xn).gt.cells(i+1,j,k,xn)) then
					Ax=-1.d0*xArea
				else if (cells(i,j,k,xn).eq.cells(i+1,j,k,xn)) then
					Ax=0.d0
				endif

				Avec(1:3)=(/Ax,Ay,Az/)

!		Setup normal vector, norm = Avec/|Avec|:
				magA = DSQRT(Avec(1)**2+Avec(2)**2+Avec(3)**2)
				norm(1:3)=(/Avec(1)/magA,Avec(2)/magA,Avec(3)/magA/)


!			Call cell select to determine velocities, pressures, and densities
				call cell_select  (xyz,i,j,k,norm,d1,u1x,u1vec,p1,d2,u2x,u2vec,p2)



				If (orderBool.eqv..FALSE.) then
!		First order evaluation has been selected:
!				Since solving in x-direction:
					uL=u1x
					uR=u2x
				
					uLvec = u1vec
					uRvec = u2vec

					dL=d1
					pL=p1
					dR=d2
					pR=p2

				else if (orderBool.eqv..TRUE.) then
!		Second order evalutaion has been selected:			

!				Extrapolate density, velocity, and pressure from gradients:
					call extrapolate (xyz,i,j,k,d1,u1vec,p1,d2,u2vec,p2,&
							& dL,uLvec,pL,dR,uRvec,pR,time)

!				Since solving in x-direction:
					uL=uLvec(1)
					uR=uRvec(1)
				endif

!			Declare location indicator for error subroutine:
			!Checks for zero or negative densities and pressures from Cell Select to prevent NaNs:
!			error_sub(IO_error,error_loc,NcXrdum,NcYrdum,NcZrdum,gamdum,m,E,V,p1,p2,d1,d2,i,j,k,time)
				error_loc="CRCS"
				call error_sub(0,error_loc,NcXr,NcYr,NcZr,gam,1.d0,1.d0,1.d0,pL,pR,dL,dR,i,j,k,time)


!###		VVV Solve Riemann Problem VVV			######


!			Solve with Riemann solver:
				call riemann(dL, uL, pL, dR, uR, pR, gg, ds,usx, ps, uflag, x_over_t, i, j, k, xyz, time)

!			Declare location indicator for error subroutine:
			!Checks for zero or negative densities and pressures from Riemann Solver to prevent NaNs:
!			error_sub(IO_error,error_loc,NcXrdum,NcYrdum,NcZrdum,gamdum,m,E,V,p1,p2,d1,d2,i,j,k,time)
				error_loc="CRRm"
				call error_sub(0,error_loc,NcXr,NcYr,NcZr,gam,1.d0,1.d0,1.d0,ps,1.d0,ds,1.d0,i,j,k,time)


!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)	

!			Update full velocity solution using uflag:
!				+1=>Right State
!				-1=>Left State
				If(uflag==1) then
					usvec=uRvec+(usx-uR)*norm
				else if (uflag==-1) then
					usvec=uLvec+(usx-uL)*norm
				else if (uflag==0) then
!					Vacuum!:
					usvec(1:3) = (/0.d0,0.d0,0.d0/)
				endif

!			Call calculate flux and update cells:
				call flux_and_update(time,xyz,i,j,k,Avec,usvec,ds,ps)

!------------------------------------------------------------------------------
!			Select y direction:				
				xyz='y'

!		Setup Unit Area Vector A-Vec. Only y-direction non-zero for y-axis flux:
				Ax=0.d0
				Az=0.d0
		!		Check left to right y positions:
				If (cells(i,j,k,yn).lt.cells(i,j+1,k,yn)) then
					Ay=1.d0*yArea
				else if (cells(i,j,k,yn).gt.cells(i,j+1,k,yn)) then
					Ay=-1.d0*yArea
				else if (cells(i,j,k,yn).eq.cells(i,j+1,k,yn)) then
					Ay=0.d0
				endif

				Avec(1:3)=(/Ax,Ay,Az/)

!		Setup normal vector, norm = Avec/|Avec|:
				magA = DSQRT(Avec(1)**2+Avec(2)**2+Avec(3)**2)
				norm(1:3)=(/Avec(1)/magA,Avec(2)/magA,Avec(3)/magA/)

!			Call cell select to determine velocities, pressures, and densities
				call cell_select  (xyz,i,j,k,norm,d1,u1y,u1vec,p1,d2,u2y,u2vec,p2)



				If (orderBool.eqv..FALSE.) then
!		First order evaluation has been selected:
!				Since solving in y-direction:
					uL=u1y
					uR=u2y
				
					uLvec = u1vec
					uRvec = u2vec

					dL=d1
					pL=p1
					dR=d2
					pR=p2

				else if (orderBool.eqv..TRUE.) then
!		Second order evalutaion has been selected:			
!				Extrapolate density, velocity, and pressure from gradients:
					call extrapolate (xyz,i,j,k,d1,u1vec,p1,d2,u2vec,p2,&
							& dL,uLvec,pL,dR,uRvec,pR,time)

!				Since solving in y-direction:
					uL=uLvec(2)
					uR=uRvec(2)
				endif



!			Declare location indicator for error subroutine:
			!Checks for zero or negative densities and pressures from Cell Select to prevent NaNs:
!			error_sub(IO_error,error_loc,NcXrdum,NcYrdum,NcZrdum,gamdum,m,E,V,p1,p2,d1,d2,i,j,k,time)
				error_loc="CRCS"
				call error_sub(0,error_loc,NcXr,NcYr,NcZr,gam,1.d0,1.d0,1.d0,pL,pR,dL,dR,i,j,k,time)


!###		VVV Solve Riemann Problem VVV			######


!			Solve with Riemann solver:
				call riemann(dL, uL, pL, dR, uR, pR, gg, ds, usy, ps, uflag, x_over_t, i, j, k, xyz, time)

!			Declare location indicator for error subroutine:
			!Checks for zero or negative densities and pressures from Riemann Solver to prevent NaNs:
!			error_sub(IO_error,error_loc,NcXrdum,NcYrdum,NcZrdum,gamdum,m,E,V,p1,p2,d1,d2,i,j,k,time)
				error_loc="CRRm"
				call error_sub(0,error_loc,NcXr,NcYr,NcZr,gam,1.d0,1.d0,1.d0,ps,1.d0,ds,1.d0,i,j,k,time)


!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)	

!			Update full velocity solution using uflag:
!				+1=>Right State
!				-1=>Left State
				If(uflag==1) then
					usvec=uRvec+(usy-uR)*norm
				else if (uflag==-1) then
					usvec=uLvec+(usy-uL)*norm
				else if (uflag==0) then
!					Vacuum!:
					usvec(1:3) = (/0.d0,0.d0,0.d0/)
				endif

!			Call calculate flux and update cells:
				call flux_and_update(time,xyz,i,j,k,Avec,usvec,ds,ps)

!!------------------------------------------------------------------------------
!			Select z direction:
				xyz="z"
				

!		Setup Unit Area Vector A-Vec. Only z-direction non-zero for z-axis flux:
				Ax=0.d0
				Ay=0.d0
		!		Check left to right y positions:
				If (cells(i,j,k,zn).lt.cells(i,j,k+1,zn)) then
					Az=1.d0*zArea
				else if (cells(i,j,k,zn).gt.cells(i,j,k+1,zn)) then
					Az=-1.d0*zArea
				else if (cells(i,j,k,zn).eq.cells(i,j,k+1,zn)) then
					Az=0.d0
				endif

				Avec(1:3)=(/Ax,Ay,Az/)

!		Setup normal vector, norm = Avec/|Avec|:
				magA = DSQRT(Avec(1)**2+Avec(2)**2+Avec(3)**2)
				norm(1:3)=(/Avec(1)/magA,Avec(2)/magA,Avec(3)/magA/)

!			Call cell select to determine velocities, pressures, and densities
				call cell_select  (xyz,i,j,k,norm,d1,u1z,u1vec,p1,d2,u2z,u2vec,p2)



				If (orderBool.eqv..FALSE.) then
!		First order evaluation has been selected:
!				Since solving in z-direction:
					uL=u1z
					uR=u2z
				
					uLvec = u1vec
					uRvec = u2vec

					dL=d1
					pL=p1
					dR=d2
					pR=p2

				else if (orderBool.eqv..TRUE.) then
!		Second order evalutaion has been selected:			
!				Extrapolate density, velocity, and pressure from gradients:
					call extrapolate (xyz,i,j,k,d1,u1vec,p1,d2,u2vec,p2,&
							& dL,uLvec,pL,dR,uRvec,pR,time)

!				Since solving in z-direction:
					uL=uLvec(3)
					uR=uRvec(3)
				endif



!			Declare location indicator for error subroutine:
			!Checks for zero or negative densities and pressures from Cell Select to prevent NaNs:
!			error_sub(IO_error,error_loc,NcXrdum,NcYrdum,NcZrdum,gamdum,m,E,V,p1,p2,d1,d2,i,j,k,time)
				error_loc="CRCS"
				call error_sub(0,error_loc,NcXr,NcYr,NcZr,gam,1.d0,1.d0,1.d0,pL,pR,dL,dR,i,j,k,time)


!###		VVV Solve Riemann Problem VVV			######


!			Solve with Riemann solver:
				call riemann(dL, uL, pL, dR, uR, pR, gg, ds, usz, ps, uflag, x_over_t, i, j, k, xyz, time)

!			Declare location indicator for error subroutine:
			!Checks for zero or negative densities and pressures from Riemann Solver to prevent NaNs:
!			error_sub(IO_error,error_loc,NcXrdum,NcYrdum,NcZrdum,gamdum,m,E,V,p1,p2,d1,d2,i,j,k,time)
				error_loc="CRRm"
				call error_sub(0,error_loc,NcXr,NcYr,NcZr,gam,1.d0,1.d0,1.d0,ps,1.d0,ds,1.d0,i,j,k,time)


!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)	

!			Update full velocity solution using uflag:
!				+1=>Right State
!				-1=>Left State
				If(uflag==1) then
					usvec=uRvec+(usz-uR)*norm
				else if (uflag==-1) then
					usvec=uLvec+(usz-uL)*norm
				else if (uflag==0) then
!					Vacuum!:
					usvec(1:3) = (/0.d0,0.d0,0.d0/)
				endif

!			Call calculate flux and update cells:
				call flux_and_update(time,xyz,i,j,k,Avec,usvec,ds,ps)

				i=i+1
			enddo
			j=j+1
		enddo
		k=k+1
	enddo

end subroutine cells_to_riemann
!-------------------------------------------------------------------------------
! Subroutine to select cell and relevant pressure, density, and velocity components
! dependent on which component is being selected, and what uflag is from previous
! component Riemann solve.
!-------------------------------------------------------------------------------
subroutine cell_select (xyz,i,j,k,norm,d1,u1cmpnt,u1vec,p1,d2,u2cmpnt,u2vec,p2)
	use constants
	implicit none

		Character(len=1),intent(in) :: xyz				!Character passed in from cells to Riemann that tells this subroutine whether we are looking at x, y, or z direction.
		Integer, intent(in) :: i,j,k					!Integers from Cells to that identify which cells the calculation needs data for
		Real*8, dimension(3),intent(in)  :: norm
		Real*8, intent (out) :: d1,p1,d2,p2
		Real*8, intent (out) :: u1cmpnt,u2cmpnt				!Component of velocity produce by dot product of norm and u1vec, u2vec
		Real*8, dimension(3), intent (out) :: u1vec,u2vec				!vector form of left and right cell velocities


!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)	

	If(xyz=="x") then
!	Here we set all left cell to be the cell we are looking at, and the right cell to be the one 
!	to the right.
		d1= cells(i,j,k,rhon)
		u1vec=(/cells(i,j,k,uxn),cells(i,j,k,uyn),cells(i,j,k,uzn)/)
		u1cmpnt=norm(1)*u1vec(1)+norm(2)*u1vec(2)+norm(3)*u1vec(3)
		p1= cells(i,j,k,pn)

		d2= cells(i+1,j,k,rhon)
		u2vec=(/cells(i+1,j,k,uxn),cells(i+1,j,k,uyn),cells(i+1,j,k,uzn)/)
		u2cmpnt=norm(1)*u2vec(1)+norm(2)*u2vec(2)+norm(3)*u2vec(3)
		p2= cells(i+1,j,k,pn)


	else if (xyz=="y") then
!	Here we set all left cell to be the cell we are looking at, and the right cell to be the one
!	in the "right" y-direction
		d1= cells(i,j,k,rhon)
		u1vec=(/cells(i,j,k,uxn),cells(i,j,k,uyn),cells(i,j,k,uzn)/)
		u1cmpnt=norm(1)*u1vec(1)+norm(2)*u1vec(2)+norm(3)*u1vec(3)
		p1= cells(i,j,k,pn)

		d2= cells(i,j+1,k,rhon)
		u2vec=(/cells(i,j+1,k,uxn),cells(i,j+1,k,uyn),cells(i,j+1,k,uzn)/)
		u2cmpnt=norm(1)*u2vec(1)+norm(2)*u2vec(2)+norm(3)*u2vec(3)
		p2= cells(i,j+1,k,pn)			

	else if (xyz=="z") then
!	Here we set all left cell to be the cell we are looking at, and the right cell to be the one
!	in the "right" y-direction
		d1= cells(i,j,k,rhon)
		u1vec=(/cells(i,j,k,uxn),cells(i,j,k,uyn),cells(i,j,k,uzn)/)
		u1cmpnt=norm(1)*u1vec(1)+norm(2)*u1vec(2)+norm(3)*u1vec(3)
		p1= cells(i,j,k,pn)

		d2= cells(i,j,k+1,rhon)
		u2vec=(/cells(i,j,k+1,uxn),cells(i,j,k+1,uyn),cells(i,j,k+1,uzn)/)
		u2cmpnt=norm(1)*u2vec(1)+norm(2)*u2vec(2)+norm(3)*u2vec(3)
		p2= cells(i,j,k+1,pn)		

	endif

end subroutine cell_select

!-------------------------------------------------------------------------------
! Subroutine to calculate flux exchange, and update cells. This uses the solutions
! calculated by the Riemann solver in cells_to_Riemann.
!-------------------------------------------------------------------------------
subroutine flux_and_update(time,xyz,i,j,k,Avec,usvec,ds,ps)
	use constants
	implicit none

		Real*8, intent(in) :: time
		Character(len=1),intent(in) :: xyz
		Real*8, dimension(3),intent(in) :: Avec, usvec
		Real*8, intent(in) :: ds,ps
		Integer, intent(in) :: i,j,k
		Real*8, dimension(3) :: Fmvec,FEvec					!Vectors to hold mass and energy flux
		Real*8, dimension(3,3) :: Fqtens					!Tensor of momentum flux
		Real*8, dimension(3,3) :: Tutens					!Tensor of Uvec Uvec^T. Acts as intermediary in calculation
		Real*8, dimension (3,3) :: unitMat					!3x3 unit matrix		
		Real*8, dimension (3) :: Fqx,Fqy,Fqz					!x,y,z components of momentum flux tensor
		Real*8 :: AFqx,AFqy,AFqz,AFm,AFE					!A dot Fq x,y,z; A dot Fm; A dot FE
		Integer :: h,l,r,c
		Real*8 :: usvec_sqrd,Fq							!The square of usVector, i.e. the dot product of us with itself
		Real*8 :: m1,m2,E1,E2,q1x,q2x,q1y,q2y,q1z,q2z
		Integer :: iL,iR,jL,jR,kL,kR						!Left and right cell indices which are varied dependent on axis being analysed.
		Character(len=4) :: error_loc						!Error Location character string· Identifies error location for error subroutine

!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)	

	usvec_sqrd = (usvec(1)*usvec(1)+usvec(2)*usvec(2)+usvec(3)*usvec(3))

!	Calculate mass and energy flux for each component in turn:
	h=1
	Do while(h<=3)
		Fmvec(h) = ds*usvec(h)
		FEvec(h) = ( (gam/(gam-1.d0))*ps + 0.5d0*ds*(usvec_sqrd) )*usvec(h)

		h=h+1
	enddo

!	Calculate Tensor momentum flux. This is a vector again when the in components,
!	and becomes a scalar when the component has a dot product taken with the area unit vector A-Vec.
!	This then gives the three components of momentum transfer.

!		First, setup unit matrix:
	h=1
	l=1
	Do while (h<=3)
		l=1
		Do while (l<=3)
			If(h==l) then
				unitMat(h,l)=1.d0
			else
				unitMat(h,l)=0.d0
			endif
			l=l+1
		enddo
		h=h+1
	enddo

!	Then, setup velocity velocity tensor:
	r=1
	c=1
	Do while (c<=3)
		r=1
		Do while (r<=3)
			Tutens(r,c) = usvec(r)*usvec(c)
			r=r+1
		enddo
		c=c+1
	enddo

!	Fill total momentum tensor:
!	Fq// = rho*(u/)(u/) + p (1//)_3	
	FqTens=(ds*Tutens+ps*unitMat)

!	Take x,y, and z components of this flux tensor:
	Fqx(1:3) = FqTens(1:3,1)!(/ds*(usvec(1)*usvec(1))+ps,ds*(usvec(1)*usvec(2)),ds*(usvec(1)*usvec(3))/)!FqTens(1:3,1)
	Fqy(1:3) = FqTens(1:3,2)!(/ds*(usvec(2)*usvec(1)),ds*(usvec(2)*usvec(2))+ps,ds*(usvec(2)*usvec(3))/)!FqTens(1:3,2)
	Fqz(1:3) = FqTens(1:3,3)!(/ds*(usvec(3)*usvec(1)),ds*(usvec(3)*usvec(2)),ds*(usvec(3)*usvec(3))+ps/)!FqTens(1:3,3)

!	Take dot product between flux tensor components and A:
	AFqx = (Avec(1)*Fqx(1))+(Avec(2)*Fqx(2))+(Avec(3)*Fqx(3))
	AFqy = (Avec(1)*Fqy(1))+(Avec(2)*Fqy(2))+(Avec(3)*Fqy(3))
	AFqz = (Avec(1)*Fqz(1))+(Avec(2)*Fqz(2))+(Avec(3)*Fqz(3))

!	Take dot products of A with mass flux and energy flux vectors:
	AFm = (Avec(1)*Fmvec(1))+(Avec(2)*Fmvec(2))+(Avec(3)*Fmvec(3))
	AFE = (Avec(1)*FEvec(1))+(Avec(2)*FEvec(2))+(Avec(3)*FEvec(3))

!	Update cells with new primitive variables:
	If (xyz=="x") then
!	Select x-axis cells to vary. Cell in right on this axis is cell number +1.
	iL=i
	iR=i+1
	jL=j
	jR=j
	kL=k
	kR=k
	else if (xyz=="y") then
!	Select y-axis cells to vary. Cell in right on this axis is cell number +1.
	iL=i
	iR=i
	jL=j
	jR=j+1
	kL=k
	kR=k
	else if(xyz=="z") then
!	Select z-axis cells to vary. Cell in right on this axis is cell number +1.
	iL=i
	iR=i
	jL=j
	jR=j
	kL=k
	kR=k+1
	else
		Print*,"Axis indicator failure [@Flux and Update]! Program will Abort!"
		STOP
	endif

!	Calculate primitive variable flux transfers:
	m1=cells(iL,jL,kL,mn) - AFm*delta_t
	q1x=cells(iL,jL,kL,qxn) - AFqx*delta_t
	q1y=cells(iL,jL,kL,qyn) - AFqy*delta_t
	q1z=cells(iL,jL,kL,qzn) - AFqz*delta_t
	E1=cells(iL,jL,kL,En) - AFE*delta_t

	m2=cells(iR,jR,kR,mn) + AFm*delta_t
	q2x=cells(iR,jR,kR,qxn) + AFqx*delta_t
	q2y=cells(iR,jR,kR,qyn) + AFqy*delta_t
	q2z=cells(iR,jR,kR,qzn) + AFqz*delta_t
	E2=cells(iR,jR,kR,En) + AFE*delta_t
	
!	Update cells for given axis Riemann solve:
	Cells(iL,jL,kL,mn)=m1
	cells(iL,jL,kL,qxn)=q1x
	cells(iL,jL,kL,qyn)=q1y
	cells(iL,jL,kL,qzn)=q1z
	cells(iL,jL,kL,En)=E1

	Cells(iR,jR,kR,mn)=m2
	cells(iR,jR,kR,qxn)=q2x
	cells(iR,jR,kR,qyn)=q2y
	cells(iR,jR,kR,qzn)=q2z
	cells(iR,jR,kR,En)=E2


!	Declare location indicator for error subroutine:
!	Checks for zero or negative masses or energies of left and then right variables, to prevent NaNs:
!	error_sub(IO_error,error_loc,NcXrdum,NcYrdum,NcZrdum,gamdum,m,E,V,p1,p2,d1,d2,i,j,k,time)

! error_sub(IO_error,error_loc,NcXrdum,NcYrdum,NcZrdum,gamdum,m,E,V,p1,p2,d1,d2,i,j,k,time)
	error_loc="flxL"
	call error_sub(0,error_loc,NcXr,NcYr,NcZr,gam,m1,E1,1.d0,1.d0,1.d0,1.d0,1.d0,i,j,k,time)

	error_loc="flxR"
	call error_sub(0,error_loc,NcXr,NcYr,NcZr,gam,m2,E2,1.d0,1.d0,1.d0,1.d0,1.d0,i,j,k,time)
end subroutine flux_and_update

!-------------------------------------------------------------------------------
!Subroutine for error messages. Takes in time, cell numbers (i,j,k), 
!physical quantities to be tested and error
!location, and prints relevanmt
!-------------------------------------------------------------------------------
subroutine error_sub(IO_error,error_loc,NcXrdum,NcYrdum,NcZrdum,gamdum,m,E,V,p1,p2,d1,d2,i,j,k,time)
	use constants
	implicit none

		Integer, intent(in) :: IO_error					!Contains input-output error value. Equals zero if no-error
		Character(len=4),intent(in) :: error_loc			!Contains character string identifying error location
		Integer, intent(in) :: NcXrdum,NcYrdum,NcZrdum			!Dummy variables for grid dimensions
		Real*8, intent(in) :: gamdum					!Dummy variable for gamma - specific heat ratio constant
		Real*8, intent (in) :: m,E,V,p1,p2,d1,d2			!Physical variables to be tested
		Integer, intent(in) :: i,j,k					!Cell location
		Real*8, intent(in) :: time					!Time at which the error is being analysed

!"cons" -> Initialise Constants
!"varC" -> Variable Conversion
!"CRCS" -> Cells to Riemann Cell Select
!"CRRm" -> Cells to Riemann Riemann Solver
!"flxL" -> Flux and Update Left Variables
!"flxR" -> Flux and Update Right Variables
!"extr" -> Extrapolate Variables from Gradients

!	GOTO 10 should skip program to Error Location key and stop.
!	GOTO 20 should skip the program printing the Error location key
!	and "STOP", and instead the program should exit this subroutine and continue.


	If(error_loc.eq."cons") then
!	In initialise constants, user inputs are checked for physicality and
!	for valid entries.
		If (IO_error.ne.(0)) then
			Print*,"Invalid Entry! Program will Abort!"
			Print*, "Location:",error_loc
			GOTO 10
		else if((NcXrdum<=0).or.(NcYrdum<=0).or.(NcZrdum<=0)) then
			Print*,"Non-Physical Grid Dimension Detected (N<=0)!"
			Print*,"Nx=",NcXrdum,"Ny=",NcYrdum,"Nz=",NcZrdum
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if (time<=0.d0) then
			Print*,"Non-physical Time Entry (t<=0)!"
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if (gam.eq.1.d0) then
!		Stability criterion to prevent gamma singularities. If gamma singularities
!& 		occur, then program will stop.
			Print*,"Gamma=1.0! Singularities detected. Program will abort!"
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if (gam.gt.(5.d0/3.d0)) then
			Print*,"Gamma>5/3 is unphysical! Program will abort!"
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else
			GOTO 20
		endif

	else if (error_loc.eq."varC") then
		If (m<0.d0) then
			Print*,"Non-Physical Mass Detected! (m<=0)"
			Print*,"m=",m
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if ((E<0.d0).and.(IsoBool .eqv. .FALSE.)) then
			Print*,"Non-Physical Energy Detected! (E<0)"
			Print*,"E=",E
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if (V<=0.d0) then
			Print*,"Non-Physical Volume Detected! (V<0)"
			Print*,"V=",V
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if (p1<0.d0) then
			Print*,"Non-Physical Pressure Detected! (P<0)"
			Print*,"P=",p1
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else
			GOTO 20
		endif

	else if (error_loc.eq."CRCS") then
		If ((d1<0.d0).or.(d2<0.d0)) then		
			Print*,"Non-Physical Density Detected! (Rho<0)"
			Print*,"Rho1=",d1,"Rho2=",d2
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if ((p1<0.d0).or.(p2<0.d0)) then		
			Print*,"Non-Physical Pressure Detected! (P<0)"
			Print*,"P1=",p1,"P2=",p2
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else
			GOTO 20
		endif

	else if (error_loc.eq."CRRm") then
		If (d1<0.d0) then		
			Print*,"Non-Physical Density Detected! (Rho<0)"
			Print*,"Rho_Sol=",d1
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if (p1<0.d0) then		
			Print*,"Non-Physical Pressure Detected! (P<0)"
			Print*,"P_Sol=",p1
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else
			GOTO 20
		endif

	else if ((error_loc.eq."flxL").or.(error_loc.eq."flxR")) then
		If (m<0.d0) then
			Print*,"Non-Physical Mass Detected! (m<0)"
			Print*,"m=",m
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if ((E<0.d0).and.(IsoBool .eqv. .FALSE.)) then
			Print*,"Non-Physical Energy Detected! (E<0)"
			Print*,"E=",E
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else
			GOTO 20
		endif
	else if (error_loc.eq."extr") then
		If ((p1<0.d0).or.(p2<0.d0)) then		
			Print*,"Non-Physical Pressure Detected! (P<0)"
			Print*,"P1=",p1,"P2=",p2
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if ((d1<0.d0).or.(d2<0.d0)) then		
			Print*,"Non-Physical Density Detected! (Rho<0)"
			Print*,"Rho1=",d1,"Rho2=",d2
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else
			GOTO 20
		endif
	else
		Print*,"Location identifier corrupted! Error Check cannot continue."
		Print*,"Error location Identifier:",error_loc
		Print*,"Program will Abort!"
		GOTO 10
	endif



10 	CONTINUE

	Print*,
	Print*,"Error Location Key:"
	Print*,"cons -> Initialise Constants"
	Print*,"varC -> Variable Conversion"
	Print*,"CRCS -> Cells to Riemann Cell Select"
	Print*,"CRRm -> Cells to Riemann Riemann Solver"
	Print*,"flxL -> Flux and Update Left Variables"
	Print*,"flxR -> Flux and Update Right Variables"
	Print*,"extr -> Extrapolate Variables from Gradients"

	STOP

20	CONTINUE

end subroutine error_sub
