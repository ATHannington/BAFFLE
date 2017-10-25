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
!* @file init_conds_squareTest.f90
!*
!* @Initial Conditions subroutine for 2D Square Test, to test viscosity, diffusion
!* and general behaviour of BAFFLE code: Fortran 90 version.
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
		Real*8 :: delx,dely,delz					!Ranges in x, y and z for the maximum range limit

!	Cells format: 	cells(i,j,k,1:14) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz/)
!	uy=0.d0
	uz=0.d0


!Sets values of pressure, x-velocity, and density according to starting conditions for 2D Sod Shock
	If ((x>=0.25d0).and.(x<=0.75d0).and.(y>=0.25d0).and.(y<=0.75d0)) then
		p=2.5d0
		rho=4.0d0
	else
		p=2.5d0
		rho=1.0d0
	endif
	
	ux=0.5d0
	uy=0.5d0

end subroutine init_conds
!===============================================================================
