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
!* @file init_conds_KH.f90
!*
!* @Initial Conditions subroutine for 2D Kelvin Helmholtz instability Simulation
!* : Fortran 90 version.
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
		Real*8 :: d,A,sig_sqrd,freq_arg						!Middle layer thickness; constants A and sigma
		Real*8 :: exparg1, exparg2, uxarglow,uxarghigh			!arguments for functions in velocities


!	uy=0.d0
	uz=0.d0

!	Set constants for Kelvin-Helmholtz:
	d=0.025d0!0.025d0
	A=0.1d0!0.1d0
	sig_sqrd= 0.00125d0!0.0125d0!0.005d0/Dsqrt(2.0d0)!0.00125d0
	
	freq_arg=1.d0

!Sets values of density according to starting conditions for 2D Kelvin-Helmholtz
	If ((y>=0.d0).and.(y<0.25d0)) then
		rho=1.d0
	else if ((y>=0.25d0).and.(y<=0.75d0)) then
		rho=10.d0
	else if ((y>0.75d0).and.(y<=1.0d0)) then
		rho=1.d0
	else
		Print*,"Invalid coordinate [density]! Program will terminate!"		!Exits program if y lies outside of 1 to 0 range
		STOP	
	endif

!Set ux arguments:
	uxarglow = ((y+d-0.25d0)/(2.d0*d))
	uxarghigh = ((y+d-0.75d0)/(2.d0*d))

!Set X-Velocity for K-H:
	If ((y>=0.d0).and.(y<=(0.25d0-d)))then
		ux=-0.5d0
	else if ((y>(0.25d0-d)).and.(y<(0.25d0+d))) then
		ux=-0.5d0+uxarglow
	else if ((y>=(0.25d0+d)).and.(y<=(0.75d0-d))) then
		ux=0.5d0
	else if ((y>(0.75d0-d)).and.(y<(0.75d0+d))) then
		ux=0.5d0-uxarghigh
	else if ((y>=0.75d0+d).and.(y<=1.d0)) then
		ux=-0.5d0
	else

!	If ((y>=0.d0).and.(y<0.25d0)) then
!		ux=-0.5d0
!	else if ((y>=0.25d0).and.(y<=0.75d0)) then
!		ux=0.5d0
!	else if ((y>0.75d0).and.(y<=1.d0)) then
!		ux=-0.5d0
!	else
		Print*,"Invalid coordinate [Vx]! Program will terminate!"		!Exits program if y lies outside of 1 to 0 range
		STOP	
	endif

!Set y-velocity:

	exparg1= -1.d0*(((y-0.25d0)**2)/(2.d0*(sig_sqrd)))
	exparg2= -1.d0*(((y-0.75d0)**2)/(2.d0*(sig_sqrd)))

	uy=(A*DSIN(freq_arg*4.d0*Pi*x))*(DEXP(exparg1)+DEXP(exparg2))

	p=2.5d0
end subroutine init_conds
!===============================================================================
