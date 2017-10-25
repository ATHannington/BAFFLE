!********************************************************************************
!* This file is part of:
!* BAFFLE (https://github.com/ATHannington/BAFFLE)
!* BAFFLE: BAsic Fortran FLuids Engine
!*
!* BAFFLE:
!* Copyright (C) 2017 Andrew Hannington (ath4@st-andrews.ac.uk)
!*
!* This software is a derivative of work by Bert Vandenbroucke
!* (bert.vandenbroucke@gmail.com)
!* Find more work by Bert Vandenbroucke at: (https://github.com/bwvdnbro)
!*
!* BAFFLE is a free software: you can redistribute it and/or modify it under the 
!* terms of the GNU Affero General Public License
!* as published by the Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* BAFFLE is distributed in the hope that it will 
!* be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!* GNU Affero General Public License for more details.
!*
!* You should have received a copy of the GNU Affero General Public License
!* along with BAFFLE. If not, see
!* <http://www.gnu.org/licenses/>.
!********************************************************************************

!********************************************************************************
!* @file hydrocode_2D.f90
!*
!* @MAIN Program file for 2D BAFFLE: Fortran 90 version
!*
!* @author Andrew Hannington (ath4@st-andrews.ac.uk)
!********************************************************************************

! Andrew Hannington								
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

!	Initialise cells array:
	call cell_setup

!	Apply boundary conditions:
	call bounds

!	Performs time integral using primitive variable conversion and riemann
!	solver, from time t=0 to time t=0.2
	call time_int

!	Write Data at end only if user has selected not to print FRAMES:
	If (writeFramesBool==0) then
		call write_data(0)
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
				  & CsQ
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

	If(qBool==1) then
	!	The following section allows the user to select second or first order evaluation:
		Print*,"Select first (f) or second (s) order evaluation:"
		Read*, orderQ

	!		Set global boolean to contain this choice:
		If ((orderQ.eq."f").or.(orderQ.eq."F")) then
			orderBool=0
		else if ((orderQ.eq."s").or.(orderQ.eq."S")) then
			orderBool=1
		else
			Print*,"Order Select Boolean Not Selected! Default Set:",orderBooldefault
			Print*,"0 => First Order	1 => Second Order"
			orderBool=orderBoolDefault
		endif

	else if (qBool==0) then
		orderBool=orderBooldefault
	endif



	If(qBool==1) then
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
			Print*,NcXrdefault,NcYrdefault,NcZrdefault
			NcXr=NcXrdefault
			NcYr=NcYrdefault
			NcZr=1.d0!NcZrdefault
		endif

	else if (qBool==0) then
		NcXr=NcXrdefault
		NcYr=NcYrdefault
		NcZr=1.d0!NcZrdefault
	endif


	NcX=NcXr+Ncellsg
	NcY=NcYr+Ncellsg
	NcZ=NcZr!+Ncellsg

	Ncells=NcX*NcY*NcZ							!Number of cells is equal to number of real cells plus number of ghost boundary cells
	allocate (cells(NcX,NcY,NcZ,29))					!Allocate the number of cells the space is divided into (Ncells) to the number of 
										!& rows to be used (cells), and the number of variables to be stored (23) to the number of columns.

!	Sets ranges of x, y and z coordinates				
	xrangeMax=1.d0
	xrangeMin=0.d0

	yrangeMax=1.d0
	yrangeMin=0.d0

	zrangeMax=0.d0
	zrangeMin=0.d0


	If(qBool==1) then
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

	else if(qBool==0) then
		gam=gamDefault
	endif



!	Suggests end-time for supersonic wind and blobs calculation:

x_Dist =0.75d0
t_Suggested = (x_Dist*(xRangemax-xRangeMin))/(1.5d0*(DSQRT(gam*(0.1d0/0.5d0))))

Print*,"Suggested End Time (Time for winds to travel ",x_Dist," Lengths of System):"
Print*,t_Suggested




	If(qBool==1) then
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

	else if(qBool==0) then

		t=tdefault

	endif

	If(qBool==1) then

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

	else if(qBool==0) then
	
		C_CFL=CFLDefault

	endif


!	Question to select isothermal or adiabatic settings

	If(qBool==1) then

		Print*,"Do you want Isothermal (i) or Adiabatic (a) settings?"
		Read*, isoQ
		If ((isoQ=="i").or.(isoQ=="I")) then
			IsoBool = 1
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
			IsoBool=0
			Print*,"ADIABATIC Selected!"
		else
			Print*,"Default set:",IsoBoolDefault
			Print*," 1=> Isothermal Setup ON	0 => Adiabatic Setup ON"
			IsoBool=IsoBoolDefault

			If(IsoBool==1) then
				Print*,"ISOTHERMAL: Default Global Cs set:",CsIsoDefault
				CsIso = CsIsoDefault
			endif
		endif

	else if(qBool==0) then
	
		IsoBool=IsoBoolDefault
		CsIso = CsIsoDefault

	endif

!	Set Initial timestep:
	delta_t=delt_default

! 	Set dummy sound speed to prevent machine precision being used in error:
	Cs=1.d0


!	Cell volume:
	CellVol = ((xrangeMax-xrangeMin)/NcXr)*((yrangeMax-yrangeMin)/NcYr)!*((zrangeMax-zrangeMin)/NcZr)

!	Set surface areas of faces
	xArea=((yrangeMax-yrangeMin)/NcYr)!*((zrangeMax-zrangeMin)/NcZr)
	yArea=((xrangeMax-xrangeMin)/NcYr)!*((zrangeMax-zrangeMin)/NcZr)
	zArea=1.d0!((xrangeMax-xrangeMin)/NcXr)*((yrangeMax-yrangeMin)/NcYr)

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
	Do while (k<=1)
		j=1
		Do while (j<=NcY)
			i=1
			Do while (i<=NcX)

				If (i==1) then
					Midx=(delx/(NcX))*(-0.5d0)
					Midy=(dely/(NcY))*((j-1)-0.5d0)
				Midz=0.d0!(delz/(NcZ))*((k-1)-0.5d0)!Sets the middle of the cell to be at the middle of the kth fraction of the z
				else if (j==1) then
					Midy=(dely/(NcY))*(-0.5d0)
					Midx=(delx/(NcX))*((i-1)-0.5d0)!Sets the middle of the cell to be at the middle of the ith fraction of the x range
				Midz=0.d0!(delz/(NcZ))*((k-1)-0.5d0)!Sets the middle of the cell to be at the middle of the kth fraction of the z
				else if (k==1) then
					Midx=(delx/(NcX))*((i-1)-0.5d0)!Sets the middle of the cell to be at the middle of the ith fraction of the x range
					Midy=(dely/(NcY))*((j-1)-0.5d0)!Sets the middle of the cell to be at the middle of the jth fraction of the y range
				Midz=0.d0!(delz/(NcZ))*(-0.5d0)

				else
					Midx=(delx/(NcX))*((i-1)-0.5d0)!Sets the middle of the cell to be at the middle of the ith fraction of the x range
					Midy=(dely/(NcY))*((j-1)-0.5d0)!Sets the middle of the cell to be at the middle of the jth fraction of the y range
				Midz=0.d0!(delz/(NcZ))*((k-1)-0.5d0)!Sets the middle of the cell to be at the middle of the kth fraction of the z range			
				endif

!			Sets values of zero for boundary cells, as otherwise they will be left with undefined values due to init_conds not applying
!			outside of the real cell range:			
				If ((i==1).or.(i==NcX).or.(j==1).or.(j==NcY)) then!.or.(k==1).or.(k==NcZ)
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

!-------------------------------------------------------------------------------
!Compact subroutine for time integration
!-------------------------------------------------------------------------------

subroutine time_int
	use constants								!Tells program to use global constants
	implicit none

		real*8 :: time							!Counter to track the time progression. Is passed to subroutines for error message purposes.	
		Real*8 :: twopercent,timedum_prcnt				!value of two percent of final time, and dummy to track value up to two percent
		Integer :: frames					!Number of frames to be written						!Run time length of movie
		Real*8 :: del_t_update						!time at which write is performed
		Real*8 :: timedum_write						!dummy percentage complete for write loop
		Real*8 :: percent
		Integer :: i,j,k
		Integer :: frame_num						!The frame number that is to be written
		character(len=1) :: writeQ
		Real*8 :: Etot	
		Integer :: vcBool	

!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)	

!	Total number of frames to be written:
	Frames = 100

	Print*,
	Print*,"Maximum Number of Data Outputs (FRAMES):", Frames

!	Set time at which data is written:
	del_t_update = t/dble(Frames)

!	Updated time 
	time = 0.d0
	timedum_prcnt=0.d0
	timedum_write=0.d0		

!	Calculate two percent of final time:
	twopercent = t*0.02d0

!	Initialise Frame Numbers:
	frame_num = 0

	If(qBool==1) then
		Print*,"Do you want to write multiple outputs of snapshots in time (FRAMES)? (y/n)"
		Read*, writeQ

		If ((writeQ=="y").or.(writeQ=="Y")) then
			writeFramesBool=1
		else if ((writeQ=="n").or.(writeQ=="N")) then
			writeFramesBool=0
		else
			Print*,"Invalid Entry! Write Frames Default Setting used:",wFramesDefault
			Print*,"0 => No Frames Written	1 => Frames Will Be Written"
			writeFramesBool=wFramesDefault
		endif

	else if(qBool==0) then

		writeFramesBool=wFramesDefault

	endif


!	Initialise boolean for variable conversion
	vcBool=0


	If (writeFramesBool==1) then
	!	Write initial conditions
		call var_conv(time,Etot)
!Print*,"t,del_t,Etot:",time,delta_t,Etot
		call write_data(frame_num)

!	Set variable conversion to false such that variable conversion is not performed twice in first run
		vcBool=1
	endif


	Do while (time<=t)

		If (vcBool==0) then
!			Convert conserved to primitive variables:
			call var_conv(time,Etot)
		endif

!Print*,"t,del_t,Etot:",time,delta_t,Etot

!		Pass cells through Riemann Solver:
		call cells_to_riemann(time)

!		Update time counter:
		time = time + delta_t

!		Calculate percentage complete
		percent = (time/t)*100.d0

!		Update time dum tracker
		timedum_prcnt = timedum_prcnt + delta_t



		If(timedum_prcnt>=twopercent) then

			Print*,
			Print'(A21,F4.0,A1)',"Percentage complete= ",percent,"%"
		!Reset percentage tracking counter to print further 2%'s
			timedum_prcnt=0.d0

			Print*,"t=",time
		endif


		If (writeFramesBool==1) then
	!		Update Write time tracker:
			timedum_write = timedum_write + delta_t

			If (timedum_write>=del_t_update) then
	!		Increase the frame number by 1 for this new frame:
				frame_num = frame_num + 1

	!		Write the data for given frame:
				call write_data(frame_num)		

	!		Reset WRITE time tracker:
				timedum_write = 0.d0

			endif
		endif


!	Set variable conversion bool back to true such that all following timesteps perform variable conversion:
		vcBool=0

	enddo


	Print*,
	Print*,"***TOTAL*** Number of Data Outputs (FRAMES):", (Frame_num+1)

end subroutine time_int

!-------------------------------------------------------------------------------
!Subroutine to write data to output files
!-------------------------------------------------------------------------------
subroutine write_data(frame_num)
	use constants								!Tells program to use global constants
	implicit none

		integer :: i,j,k
		character(len=20) :: filedata
		character(len=20) :: format_string1,format_string2,format_string3		!Format descriptors 
		character(len=100) :: filename,tmp1,tmp2						!Filename strings: final and temporary
		character(len=43) :: path							!File path for writing file
		character(len=20) :: Left,Right
		Integer,intent(in):: frame_num
		Integer :: LengthAdjustBool

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

!	Switches adjustable file string length on or off:
!	0 => OFF		1=> ON
	LengthAdjustBool = 0


!	Setup left and right of filename:
	Left="Output_"
	Right=".dat"

!	Set path to file output directory:
	path = "/home/alpha/ath4/summer2017/HydroCode/DISO/"

!	Define format strings for later filename creation:
	If (LengthAdjustBool == 0) then 
!	If adjustable filename length is off, then a preset filename size is used: 
			format_string1 = '(A7,I5.5)'
			format_string2 = '(A12,A4)'
			format_string3 = '(A43,A16)'
	else if (LengthAdjustBool==1) then
!	If adjustable filename length is on, then the zeros prefixing the frame_number ajust
!	to the number of frames being used. 
!	I.E. <10 frames -> 01,02,03 etc.	>=10 and <100 -> 001,002,003, etc.
		If (frame_num<10) then
			format_string1 = '(A7,I2.2)'
			format_string2 = '(A9,A4)'
			format_string3 = '(A43,A13)'
		else if ((frame_num>=10).and.(frame_num<100)) then
			format_string1 = '(A7,I3.3)'
			format_string2 = '(A10,A4)'
			format_string3 = '(A43,A14)'
		else if ((frame_num>=100).and.(frame_num<1000)) then
			format_string1 = '(A7,I4.4)'
			format_string2 = '(A11,A4)'
			format_string3 = '(A43,A15)'
		else if ((frame_num>=1000).and.(frame_num<10000)) then
			format_string1 = '(A7,I5.5)'
			format_string2 = '(A12,A4)'
			format_string3 = '(A43,A16)'
		else
			Print*, "LENGTH FAILURE [@ Filename Length Adjust]! PROGRAM WILL ABORT!"	
			STOP
		endif
	else
		Print*,"LengthAdjustBool Failure [@ Filename Length Adjust Switch]! Program will abort!"
		STOP
	endif

	Write(tmp1,format_string1) Left,frame_num
	write(tmp2,format_string2) tmp1,Right
	Write(filename,format_string3) path, tmp2
	
	Open(1,file=trim(filename))

!	Insert format header to each data file:
	write(1,*) (/"   x  ","   y  ","   z  ","   m  ","  qx  ","  qy  ","  qz  ",&
			"   E  "," rho  ","  ux  ","  uy  ","  uz  ", &
			"   p  ","   V  "/)

!	i=2 so data written skips ghost cell. Similarly for ncells-1 in do loop
i=2
j=2
k=1!2

	Print*,"Writing Frame:",frame_num


	Do while(k<=1)!NcX-1)
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

	Print*,"Finished Writing Frame:",frame_num

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
		Real*8 :: Ediff,E1,E2						!Difference in internal energies

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

	Do while (k<=1)!NcZ)
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

!				rho = m/V:
				rho=m/v

!				u=q/m:
				ux = qx/m
				uy = qy/m
				uz = qz/m

		!		p = (gamma - 1)*(E/V - 1/2 * rho * u^2)
				p = (gam-1.d0)*( (E/V)-( 0.5d0*rho*((ux**2)+(uy**2)+(uz**2)) ) )

				E1 = E - (0.5d0*V*rho*((ux**2)+(uy**2)+(uz**2)))

				If (IsoBool==1) then
					p = isotherm_p(rho)

!				Adjust Total energy of cell by incorporating Isothermal Loss term.
!				The energy lost is the difference in internal energies of the adiabatic
!				and Isothermal case.
					E2 = (p/(gam-1.d0))*V
					Ediff = E1-E2
					E = E - Ediff
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
	If (orderBool.eq.1) then
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

	delmin=minval((/delx,dely/))!,delz/))

!	Calculate sound speed for given cell:
	cs = DSQRT(gam*(p/rho))	

	umag=DSQRT(ux**2+uy**2+uz**2)

!	Calculate maximum time step for given cell
	deltmax = C_CFL*(delmin/(umag+Cs))

!	For first cell only, set new minimum timestep:
	If ((i==1).and.(j==1).and.(k==1)) then
		delta_t=deltmax
!	If newly calucated timestep is less than previous smallest timestep, set timestep
!	to new, smaller value.
	else if (deltmax.lt.delta_t) then
		delta_t=deltmax
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

	Do while (k<=1)!NcZ-1)
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



				If (orderBool.eq.0) then
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

				else if (orderBool.eq.1) then
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
				call riemann(dL, uL, pL, dR, uR, pR, gg, ds, usx, ps, uflag, x_over_t)

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
					Print*,"Uflag=0.d0! Vacuum detected! Program will abort!"
					STOP
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



				If (orderBool.eq.0) then
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

				else if (orderBool.eq.1) then
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
				call riemann(dL, uL, pL, dR, uR, pR, gg, ds, usy, ps, uflag, x_over_t)

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
					Print*,"Uflag=0.d0! Vacuum detected! Program will abort!"
					STOP
				endif

!			Call calculate flux and update cells:
				call flux_and_update(time,xyz,i,j,k,Avec,usvec,ds,ps)

!!------------------------------------------------------------------------------
!!			Select z direction:
!				xyz="z"
!				
!
!!		Setup Unit Area Vector A-Vec. Only z-direction non-zero for z-axis flux:
!				Ax=0.d0
!				Ay=0.d0
!		!		Check left to right y positions:
!				If (cells(i,j,k,zn).lt.cells(i,j,k+1,zn)) then
!					Az=1.d0*zArea
!				else if (cells(i,j,k,zn).gt.cells(i,j,k+1,zn)) then
!					Az=-1.d0*zArea
!				else if (cells(i,j,k,zn).eq.cells(i,j,k+1,zn)) then
!					Az=0.d0
!				endif
!
!				Avec(1:3)=(/Ax,Ay,Az/)
!
!!		Setup normal vector, norm = Avec/|Avec|:
!				magA = DSQRT(Avec(1)**2+Avec(2)**2+Avec(3)**2)
!				norm(1:3)=(/Avec(1)/magA,Avec(2)/magA,Avec(3)/magA/)
!
!!			Call cell select to determine velocities, pressures, and densities
!				call cell_select  (xyz,i,j,k,norm,d1,u1z,u1vec,p1,d2,u2z,u2vec,p2)
!
!
!
!				If (orderBool.eq.0) then
!!		First order evaluation has been selected:
!!				Since solving in z-direction:
!					uL=u1z
!					uR=u2z
!				
!					uLvec = u1vec
!					uRvec = u2vec
!
!					dL=d1
!					pL=p1
!					dR=d2
!					pR=p2
!
!				else if (orderBool.eq.1) then
!!		Second order evalutaion has been selected:			
!!				Extrapolate density, velocity, and pressure from gradients:
!					call extrapolate (xyz,i,j,k,d1,u1vec,p1,d2,u2vec,p2,&
!							& dL,uLvec,pL,dR,uRvec,pR,time)
!
!!				Since solving in z-direction:
!					uL=uLvec(3)
!					uR=uRvec(3)
!				endif
!
!
!
!!			Declare location indicator for error subroutine:
!			!Checks for zero or negative densities and pressures from Cell Select to prevent NaNs:
!!			error_sub(IO_error,error_loc,NcXrdum,NcYrdum,NcZrdum,gamdum,m,E,V,p1,p2,d1,d2,i,j,k,time)
!				error_loc="CRCS"
!				call error_sub(0,error_loc,NcXr,NcYr,NcZr,gam,1.d0,1.d0,1.d0,pL,pR,dL,dR,i,j,k,time)
!
!
!!###		VVV Solve Riemann Problem VVV			######
!
!
!!			Solve with Riemann solver:
!				call riemann(dL, uL, pL, dR, uR, pR, gg, ds, usz, ps, uflag, x_over_t)
!
!!			Declare location indicator for error subroutine:
!			!Checks for zero or negative densities and pressures from Riemann Solver to prevent NaNs:
!!			error_sub(IO_error,error_loc,NcXrdum,NcYrdum,NcZrdum,gamdum,m,E,V,p1,p2,d1,d2,i,j,k,time)
!				error_loc="CRRm"
!				call error_sub(0,error_loc,NcXr,NcYr,NcZr,gam,1.d0,1.d0,1.d0,ps,1.d0,ds,1.d0,i,j,k,time)
!
!
!!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!!						& grad_d_x,grad_d_y,grad_d_z, &
!!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!!						& grad_p_x,grad_p_y,grad_p_z /)	
!
!!			Update full velocity solution using uflag:
!!				+1=>Right State
!!				-1=>Left State
!				If(uflag==1) then
!					usvec=uRvec+(usz-uR)*norm
!				else if (uflag==-1) then
!					usvec=uLvec+(usz-uL)*norm
!				else if (uflag==0) then
!					Print*,"Uflag=0.d0! Vacuum detected! Program will abort!"
!					STOP
!				endif
!
!!			Call calculate flux and update cells:
!				call flux_and_update(time,xyz,i,j,k,Avec,usvec,ds,ps)
!
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
!	q1z=cells(iL,jL,kL,qzn) - AFqz*delta_t
	E1=cells(iL,jL,kL,En) - AFE*delta_t

	m2=cells(iR,jR,kR,mn) + AFm*delta_t
	q2x=cells(iR,jR,kR,qxn) + AFqx*delta_t
	q2y=cells(iR,jR,kR,qyn) + AFqy*delta_t
!	q2z=cells(iR,jR,kR,qzn) + AFqz*delta_t
	E2=cells(iR,jR,kR,En) + AFE*delta_t
	
!	Update cells for given axis Riemann solve:
	Cells(iL,jL,kL,mn)=m1
	cells(iL,jL,kL,qxn)=q1x
	cells(iL,jL,kL,qyn)=q1y
!	cells(iL,jL,kL,qzn)=q1z
	cells(iL,jL,kL,En)=E1

	Cells(iR,jR,kR,mn)=m2
	cells(iR,jR,kR,qxn)=q2x
	cells(iR,jR,kR,qyn)=q2y
!	cells(iR,jR,kR,qzn)=q2z
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
		else if (time<=0) then
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
		If (m<=0) then
			Print*,"Non-Physical Mass Detected! (m<=0)"
			Print*,"m=",m
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if (E<0) then
			Print*,"Non-Physical Energy Detected! (E<0)"
			Print*,"E=",E
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if (V<=0) then
			Print*,"Non-Physical Volume Detected! (V<0)"
			Print*,"V=",V
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if (p1<=0) then
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
		If ((d1<=0.d0).or.(d2<=0.d0)) then		
			Print*,"Non-Physical Density Detected! (Rho<0)"
			Print*,"Rho1=",d1,"Rho2=",d2
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if ((p1<=0.d0).or.(p2<=0.d0)) then		
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
		If (d1<=0.d0) then		
			Print*,"Non-Physical Density Detected! (Rho<0)"
			Print*,"Rho_Sol=",d1
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if (p1<=0.d0) then		
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
		If (m<=0) then
			Print*,"Non-Physical Mass Detected! (m<=0)"
			Print*,"m=",m
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if (E<0) then
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
		If ((p1<=0.d0).or.(p2<=0.d0)) then		
			Print*,"Non-Physical Pressure Detected! (P<0)"
			Print*,"P1=",p1,"P2=",p2
			Print*,"Cell Location: i=",i,"j=",j,"k=",k
			Print*,"TimeStamp: t=",time
			Print*, "Location:",error_loc
			Print*,"Program will Abort!"
			GOTO 10
		else if ((d1<=0.d0).or.(d2<=0.d0)) then		
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
!********************************************************************************
!* This file is part of python_finite_volume_solver
!* Copyright (C) 2017 Kenneth Wood (kw25@st-andrews.ac.uk)
!*                    Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
!*
!* python_finite_volume_solver is free software: you can redistribute it and/or
!* modify it under the terms of the GN1.4999999999999999E-002U Affero General Public License as
!* published by the Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* python_finite_volume_solver is distributed in the hope that it will be useful,
!* but WITOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!* GNU Affero General Public License for more details.
!*
!* You should have received a copy of the GNU Affero General Public License
!* along with python_finite_volume_solver. If not, see
!* <http://www.gnu.org/licenses/>.
!********************************************************************************

!********************************************************************************
!* @file riemansolver.f
!*
!* @brief Stanalone Riemann solver library: Fortran 77 version.
!*
!* This Riemann solver is the original exact Riemann solver as it is presented in
!* Toro, E., Riemann Solvers and Numerical Methods for Fluid Dynamics, 3rd
!* edition (Springer, 2009), chapter 4.
!*
!* It was carefully copied (by hand) from that book to actual Fortran 77 code by
!* Kenneth Wood, and slightly modified to be used as a Riemann solver routine in
!* a 1-3D finite volume code.
!*
!* Additional changes (i1.4999999999999999E-002ncluding this documentation) were made by Bert
!* Vandenbroucke to make the Python, C++ and Fortran version of the Riemann
!* solver more homogeneous.
!*
!* @author Kenneth Wood (kw25@st-andrews.ac.uk)
!* @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
!********************************************************************************

!********************************************************************************
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
!* @param us Variable to store the resulting fluid velocity in.
!* @param ps Variable to store the resulting pressure in.
!* @param uflag Flag variable used to signal if the solution was sampled from the
!* left state (-1) or the right state (1). A vacuum state (0) is currently not
!* supported by this version of the Riemann solver.
!* @param x_over_t Point (in velocity space) where we want to sample the
!* solution. In a finite volume scheme you usually want to set this value to 0.
!********************************************************************************
      subroutine riemann(d1, u1, p1, d2, u2, p2, gg, ds,us, ps, uflag, x_over_t)

      implicit none

!     Declaration of variables:
      integer I, CELLS, uflag
      REAL*8 :: AGAMMA, G1, G2, G3, G4, G5, G6, G7, G8, &
          DL, UL, PL, CL, DR, UR, PR, CR, &
          DIAPH, DOMLEN, DS, DX, PM, MPA, PS, S, &
          TIMEOUT, UM, US, XPOS

      real*8 gg, d1, u1, p1, d2, u2, p2
      real*8 x_over_t, dxdt

      COMMON /GAMMAS/ AGAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR

      dxdt = x_over_t

!     Compute gamma related constants

      AGAMMA = gg
      DL = d1
      UL = u1
      PL = p1
      DR = d2
      UR = u2
      PR = p2

      G1 = (AGAMMA - 1.0d0)/(2.0d0*AGAMMA)
      G2 = (AGAMMA + 1.0d0)/(2.0d0*AGAMMA)
      G3 = 2.0d0*AGAMMA/(AGAMMA - 1.0d0)
      G4 = 2.0d0/(AGAMMA - 1.0d0)
      G5 = 2.0d0/(AGAMMA + 1.0d0)
      G6 = (AGAMMA - 1.0d0)/(AGAMMA + 1.0d0)
      G7 = (AGAMMA - 1.0d0)/2.0d0
      G8 = AGAMMA - 1.0d0

      MPA=1.0d0  ! Not sure about this????

!     Compute sound speeds

      CL = DSQRT(AGAMMA*PL/DL)
      CR = DSQRT(AGAMMA*PR/DR)

!     The pressure positivity condition is tested for

      IF(G4*(CL+CR).LE.(UR-UL))THEN
!     The initial data is such that vacuum is generated. Program stopped
        WRITE(6,*) 'Vacuum generated by data. Program stopped.'
        STOP
      ENDIF

!     Exact solution for pressure and velocity in star region is found

      CALL STARPU(PM, UM, MPA)
      if(pm.ne.pm) then
        print*,pm,pl,dl,ul,pr,dr,ur
        stop
      endif
      if(UM .gt.dxdt) then
        uflag=-1 ! "left state"
      else
        uflag=1 ! "right state
      endif


      CALL SAMPLE(PM, UM, dxdt, DS, US, PS)

      RETURN
      END

!********************************************************************************
!* @brief Get the pressure and velocity in the middle state.
!*
!* @param p Variable to store the resulting pressure in.
!* @param u Variable to store the resulting fluid velocity in.
!* @param mpa Mach number? Not used.
!********************************************************************************
      SUBROUTINE STARPU(P, U, MPA)

      IMPLICIT NONE

!     Purpose: to compute the solution for pressure and velocity in star region

      INTEGER I, NRITER

      REAL*8 DL, UL, PL, CL, DR, UR, PR, CR, &
         CHANGE, FL, FLD, FR, FRD, P, POLD, PSTART, &
         TOLPRE, U, UDIFF, MPA

      COMMON/STATES/ DL, UL, PL, CL, DR, UR, PR, CR
      DATA TOLPRE, NRITER/1.0d-09, 20/

!     Guessed value PSTART is computed

      CALL GUESSP(PSTART)

      POLD = PSTART
      UDIFF = UR - UL
!      WRITE(6,*)'Iteration number Change '

      DO 10 I = 1, NRITER

         CALL PREFUN(FL, FLD, POLD, DL, PL, CL)
         CALL PREFUN(FR, FRD, POLD, DR, PR, CR)
         P = POLD - (FL + FR + UDIFF)/(FLD + FRD)

         CHANGE = 2.0d0*DABS((P - POLD)/(P + POLD))
!         WRITE(6, 30)I, CHANGE
         IF(CHANGE.LE.TOLPRE)GOTO 20
         IF(P.LT.0.0)P = TOLPRE
         POLD = P

 10   CONTINUE

!      WRITE(6,*)'Divergence in Newton-Raphson iteration'

 20   CONTINUE

!     Compute velocity in star region

      U = 0.5d0*(UL + UR + FR - FL)

!      WRITE(6,*)'Pressure	Velocity'
!      WRITE(6,40)P/MPA, U

 30   FORMAT(5X, I5,15X, F12.7)
 40   FORMAT(2(F14.6, 5X))

      RETURN
      END

!********************************************************************************
!* @brief Get an initial guess for the pressure in the middle state, as a
!* starting point for the Newton-Raphson iteration.
!*
!* @param pm Variable to store the resulting pressure guess in.
!********************************************************************************
      SUBROUTINE GUESSP(PM)

!     Purpose: to provide a guess value for pressure PM in the Star Region.
!     The choice is made according to adaptive Riemann solver using the PVRS,
!     TRRS, and TSRS approximate Riemann solvers

      IMPLICIT NONE

      REAL*8 DL, UL, PL, CL, DR, UR, PR, CR, &
          AGAMMA, G1, G2, G3, G4, G5, G6, G7, G8, &
          CUP, GEL, GER, PM, PMAX, PMIN, PPV, PQ, &
          PTL, PTR, QMAX, QUSER, UM

      COMMON /GAMMAS/ AGAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR

      QUSER = 2.0d0

!     Compute guess pressure for PVRS Riemann solver

      CUP = 0.25d0*(DL + DR)*(CL + CR)
      PPV = 0.5d0*(PL + PR) + 0.5d0*(UL - UR)*CUP
      PPV = AMAX1(0.0d0, PPV)
      PMIN = AMIN1(PL, PR)
      PMAX = AMAX1(PL, PR)
      QMAX = PMAX/PMIN

      IF(QMAX.LE.QUSER.AND. &
       (PMIN.LE.PPV.AND.PPV.LE.PMAX))THEN

!        Select PVRS Riemann solver

         PM = PPV

      ELSE
         IF(PPV.LT.PMIN)THEN

!           Select Two-Rarefaction Riemann solver

            PQ = (PL/PR)**G1
            UM = (PQ*UL/CL + UR/CR + &
                 G4*(PQ - 1.0d0))/(PQ/CL + 1.0d0/CR)
            PTL = 1.0d0 + G7*(UL - UM)/CL
            PTR = 1.0d0 + G7*(UM - UR)/CR
            PM = 0.5d0*(PL*PTL**G3 + PR*PTR**G3)
         ELSE

!           Select Two=Shock Riemann solver with PVRS as estimate

            GEL = DSQRT((G5/DL)/(G6*PL + PPV))
            GER = DSQRT((G5/DR)/(G6*PR + PPV))
            PM = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER)
         ENDIF
      ENDIF

      RETURN
      END

!********************************************************************************
!* @brief Evaluate the pressure functions for the given pressure guess and the
!* given left or right state variables.
!*
!* @param f Variable to store the pressure function value in.
!* @param fd Variable to store the value of the derivative of the presure
!* function in.
!* @param p Current pressure guess.
!* @param dk Left or right state density.
!* @param pk Left or right state pressure.
!* @param ck Left or right state sound speed.
!********************************************************************************
      SUBROUTINE PREFUN(F, FD, P, DK, PK, CK)

!     Purpose: to evaluate the pressure functions FL and FR in the exact Riemann
!     solver

      IMPLICIT NONE

      REAL*8 AK, BK, CK, DK, F, FD, P, PK, PRAT, QRT, &
          AGAMMA, G1, G2, G3, G4, G5, G6, G7, G8

      COMMON /GAMMAS/ AGAMMA, G1, G2, G3, G4, G5, G6, G7, G8

      IF(P.LE.PK) THEN

!        Rarefaction wave

         PRAT = P/PK
         F = G4*CK*(PRAT**G1 - 1.0d0)
         FD = (1.0d0/(DK*CK))*PRAT**(-G2)

      ELSE

!        Shock wave
         AK = G5/DK
         BK = G6*PK
         QRT = DSQRT(AK/(BK + P))
         F = (P - PK)*QRT
         FD = (1.0d0 - 0.5d0*(P - PK)/(BK + P))*QRT

      ENDIF

      RETURN
      END

!********************************************************************************
!* @brief Sample the Riemann solution with the given middle state pressure and
!* fluid velocity at the given point in velocity space.
!*
!* @param pm Middle state pressure.
!* @param um Middle state fluid velocity.
!* @param s Sampling point in velocity space.
!* @param d Variable to store the resulting density in.
!* @param u Variable to store the resulting fluid velocity in.
!* @param p Variable to store the resulting pressure in.
!********************************************************************************
      SUBROUTINE SAMPLE(PM, UM, S, D, U, P)

!     Purpose to sample the solution throughout the wave pattern. Pressure PM
!     and velocity UM in the Star Region are known. Sampling is performed in
!     terms of the 'speed' S = X/T. Sampled values are D, U, P

!     Input variables: PM, UM, S, /GAMMAS. /STATES/
!     Output variables: D, U, P

      IMPLICIT NONE

      REAL*8 DL, UL, PL, CL, DR, UR, PR, CR, &
          AGAMMA, G1, G2, G3, G4, G5, G6, G7, G8, &
          C, CML, CMR, D, P, PM, PML, PMR, S, &
          SHL, SHR, SL, SR, STL, STR, U, UM

      COMMON /GAMMAS/ AGAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR

      IF(S.LE.UM)THEN

!     Sampling point lies to the left of the contact discontinuity

         IF(PM.LE.PL)THEN

!          Left rarefaction

            SHL = UL - CL

            IF(S.LE.SHL)THEN

!              Sampled point is left data state

               D = DL
               U = UL
               P = PL
            ELSE
               CML = CL*(PM/PL)**G1
               STL = UM - CML

               IF(S.GT.STL)THEN

!                 Sampled point is Star Left state

                  D = DL*(PM/PL)**(1.0d0/AGAMMA)
                  U = UM
                  P = PM
               ELSE

!                 Sampled point is inside left fan

                  U = G5*(CL + G7*UL + S)
                  C = G5*(CL + G7*(UL - S))
                  D = DL*(C/CL)**G4
                  P = PL*(C/CL)**G3
               ENDIF
            ENDIF
!         ENDIF
         ELSE

!            Left shock

            PML = PM/PL
            SL = UL - CL*DSQRT(G2*PML + G1)

            IF(S.LE.SL)THEN

!              Sampled point is left data state

               D = DL
               U = UL
               P = PL

             ELSE

!              Sampled point is Star Left state

               D = DL*(PML+G6)/(PML*G6 + 1.0d0)
               U = UM
               P = PM
             ENDIF
           ENDIF
          ELSE

!    Sampling point lies to the right of the contact discontinuity

         IF(PM.GT.PR)THEN

!           Right shock

            PMR = PM/PR
            SR = UR + CR*DSQRT(G2*PMR + G1)

            IF(S.GE.SR)THEN

!              Sampled point is right data state

               D = DR
               U = UR
               P = PR
            ELSE

!              Sampled point is Star Right state

               D = DR*(PMR + G6)/(PMR*G6 + 1.0d0)
               U = UM
               P = PM
            ENDIF
         ELSE


!           Right rarefaction

            SHR = UR + CR

            IF(S.GE.SHR)THEN

!              Sampled point is right data state

               D = DR
               U = UR
               P = PR
            ELSE
               CMR = CR*(PM/PR)**G1
               STR = UM + CMR

               IF(S.LE.STR)THEN

!                 Sampled point is Star Right state

                  D = DR*(PM/PR)**(1.0d0/AGAMMA)
                  U = UM
                  P = PM
               ELSE

!                 Sampled point is inside left fan

                  U=G5*(-CR + G7*UR + S)
                  C = G5*(CR - G7*(UR-S))
                  D = DR*(C/CR)**G4
                  P = PR*(C/CR)**G3
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      RETURN
END

