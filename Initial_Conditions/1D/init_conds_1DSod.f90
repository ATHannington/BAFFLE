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
!* @file init_conds_1DSod.f90
!*
!* @Initial Conditions subroutine for 1D Sod Shock Test Tube Simulation: Fortran 90 version.
!*
!* @author Andrew Hannington (ath4@st-andrews.ac.uk)
!********************************************************************************

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

		real*8,external::isotherm_p					!Isothermal Pressure
		
!#	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!#						& grad_d_x,grad_d_y,grad_d_z, &
!#						& grad_ux_x,grad_uy_x,grad_uz_x, &
!#						& grad_ux_y,grad_uy_y,grad_uz_y, &
!#						& grad_ux_z,grad_uy_z,grad_uz_z, &
!#						& grad_p_x,grad_p_y,grad_p_z /)

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
	upbound=1.d0!Sqrt(delx**2+dely**2+delz**2)

!Sets values of pressure, x-velocity, and density according to starting conditions for 2D Sod Shock
	If ((x>=lowbound).and.(x<midbound)) then
		rho=1.d0
		ux=0.d0
		uy=0.d0
		uz=0.d0
		p=1.d0
	else if ((x>=midbound).and.(x<=upbound)) then
		rho=0.125d0
		ux=0.d0
		uy=0.d0
		uz=0.d0
		p=0.1d0
	else
		Print*,"Invalid coordinate! Program will terminate!"		!Exits program if x lies outside of 1 to 0 range
		STOP	
	endif

	If (IsoBool.eqv..TRUE.) then
		p = isotherm_p(rho)
	endif

end subroutine init_conds
!===============================================================================
