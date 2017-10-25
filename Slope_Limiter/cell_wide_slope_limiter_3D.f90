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
!* @file cell_wide_slope_limiter_3D.f90
!*
!* @Seperate advanced slope limiter for slope of primitive variables calculated in 
!* get_gradients_subroutine_#.f90. 3D version : Fortran 90 version
!*
!* @author Andrew Hannington (ath4@st-andrews.ac.uk)
!********************************************************************************

!-------------------------------------------------------------------------------
!*******************************************************************************
!-------------------------------------------------------------------------------
! Implementinmg Cell-Wide slop limiter, where
!
! Grad-vec = alpha* grad_dum-vec
!
! alpha = max(0, min (1,beta*min[(phi_ngb,max - phi)/(phi_ext,max - phi), (phi - phi_ngb,min)/(phi - phi_ext,min)]))
!
!The Max(0,(~)) term ensures gradient is set to zero at turning points, as alpha!>=0
!
! phi_ngb,max = max_j[phi_j]
! phi_ngb,min = min_j[phi_j]
!
! phi_ext,max = max_j[phi+grad_dum-vec·Del_x_j-vec]
! phi_ext,min = min_j[phi+grad_dum-vec·Del_x_j-vec]
!
! phi_j is the variable (rho,ux,uy,uz,p) in neighbour (ngb) cell j
! Del_x_j -vec is the vector pointing between midpoint of the cell to the midpoint of the face between cell and neighbour
! Beta is a constant: beta=0.5 is the most conservative (stable) value allowed!
! 
!-------------------------------------------------------------------------------
!*******************************************************************************
!-------------------------------------------------------------------------------
subroutine slope_limiter(i,j,k,time)
	use constants
	implicit none

		Real*8,intent(in) :: time
		Integer,intent(in) :: i,j,k					!Incoming integers from gradients subroutine. This selects the cell who's gradient has just been calculated.
		Integer :: idum,jdum,kdum					!These dummy i,j,k integers hold the location of the neighbour cell for the specific neighbour on a given axis
		Integer :: ngb							!Identifies the number of neighbour cells.	*** 2=> 1D		4=>2d		6=>3d ***
		Real*8, allocatable,dimension(:) :: d_array,ux_array,uy_array, &
					& uz_array,p_array			!Holds the variables for the neighbour cells for the respective variable such that the min and max value can be 										determined

		Real*8 :: d_ngb_min,ux_ngb_min,uy_ngb_min,& 
			& uz_ngb_min,p_ngb_min					!Holds the minimum value for each variable IN *** NEIGHBOUR CELLS ***
		Real*8 :: d_ngb_max,ux_ngb_max,uy_ngb_max, &
			& uz_ngb_max,p_ngb_max					!Holds the maximum value for each variable IN *** NEIGHBOUR CELLS ***

		Real*8 :: d_ext_min,ux_ext_min,uy_ext_min,& 
			& uz_ext_min,p_ext_min					!Holds the minimum value for each variable *** FROM EXTRAPOLATION ***
		Real*8 :: d_ext_max,ux_ext_max,uy_ext_max, &
			& uz_ext_max,p_ext_max					!Holds the maximum value for each variable *** FROM EXTRAPOLATION ***

		Real*8,allocatable,dimension(:) :: d_ext_array,ux_ext_array, &
			& uy_ext_array, uz_ext_array,p_ext_array		!Holds the calculated value for each variable *** FROM EXTRAPOLATION *** for all neighbours, so that the min or max may 										be calculated

		Real*8 :: d_cent,ux_cent,uy_cent,& 
			& uz_cent,p_cent					!Holds the value for each variable *** IN CENTRE OF CELL ***

		Real*8, dimension(3) :: Del_x_vec				!Holds the centre to face centre vector
		Real*8 :: delx,dely,delz					!Distances for Del_x_vector between cell centre and face centre

		Real*8,dimension(3) :: g_d_vec,g_ux_vec,g_uy_vec, &		!Holds the gradient for each vector along each axis in vector form. These are the **grad_dum** vectors!!
					& g_uz_vec,g_p_vec
		Real*8 :: alpha, beta						!Alpha the scaling variable for the gradients, as defined above. Beta a constant setting how conservative the code is! 	

		Integer :: num							!Integer number that ensures all of the variables in a given direction are given the same position in arrays.

!-------------------------------------------------------------------------------
!#####				SETUP					########
!-------------------------------------------------------------------------------
!Set conservative beta:
	beta=0.5d0
!Initialise alpha as zero to begin with:
	alpha=0.d0

!Set number of neighbour cells:
!*** 2=> 1D		4=>2d		6=>3d ***
	ngb=6

!Set Dimensions of arrays that are dependent of number of neighbours:
	allocate(d_array(ngb))
	allocate(ux_array(ngb))
	allocate(uy_array(ngb))
	allocate(uz_array(ngb))
	allocate(p_array(ngb))

	allocate(d_ext_array(ngb))
	allocate(ux_ext_array(ngb))
	allocate(uy_ext_array(ngb))
	allocate(uz_ext_array(ngb))
	allocate(p_ext_array(ngb))

!Fill Gradient vector for cell being analysed:
!These are the **grad_dum** vectors!!
	g_d_vec = (/cells(i,j,k,gdxn),cells(i,j,k,gdyn),cells(i,j,k,gdzn)/)
	g_ux_vec = (/cells(i,j,k,guxxn),cells(i,j,k,guxyn),cells(i,j,k,guxzn)/)
	g_uy_vec = (/cells(i,j,k,guyxn),cells(i,j,k,guyyn),cells(i,j,k,guyzn)/)
	g_uz_vec = (/cells(i,j,k,guzxn),cells(i,j,k,guzyn),cells(i,j,k,guzzn)/)
	g_p_vec = (/cells(i,j,k,gpxn),cells(i,j,k,gpyn),cells(i,j,k,gpzn)/)

!Fill central values of the cell:
!These are the phi values:

	d_cent = cells(i,j,k,rhon)
	ux_cent = cells(i,j,k,uxn)
	uy_cent = cells(i,j,k,uyn)
	uz_cent = cells (i,j,k,uzn)
	p_cent = cells (i,j,k,pn)

!-------------------------------------------------------------------------------
!#####			X-Axis Neighbours				########
!-------------------------------------------------------------------------------

!-ve x-face direction neighbour
	idum = i-1
	jdum = j
	kdum = k

!Set position in array number:
	num=1

!Fill variable array entries for this neighbour:
! This is equivalent to entering the first value (j=1) into the phi_j, so that the min and max can be taken later.
	d_array(num) = cells(idum,jdum,kdum,rhon)
	ux_array(num) = cells(idum,jdum,kdum,uxn)
	uy_array(num) = cells(idum,jdum,kdum,uyn)
	uz_array(num) = cells(idum,jdum,kdum,uzn)
	p_array(num) = cells(idum,jdum,kdum,pn)


! *** EXTRAPOLATION ***:
!Calculate cell centre to cell face distances:
	delx = (cells(i,j,k,xn) - cells(idum,jdum,kdum,xn))/2.d0
	dely = (cells(i,j,k,yn) - cells(idum,jdum,kdum,yn))/2.d0	
	delz = (cells(i,j,k,zn) - cells(idum,jdum,kdum,zn))/2.d0	

!Fill temporary cell centre to face centre vector:
	Del_x_vec = (/delx,dely,delz/)

!Make extrapolations to cell face:
! " phi+grad_dum-vec·Del_x_j-vec "
!	d_ext_array(num) = d_cent + dot_product(g_d_vec,Del_x_vec)
!	ux_ext_array(num) = ux_cent + dot_product(g_ux_vec,Del_x_vec)
!	uy_ext_array(num) = uy_cent + dot_product(g_uy_vec,Del_x_vec)
!	uz_ext_array(num) = uz_cent + dot_product(g_uz_vec,Del_x_vec)
!	p_ext_array(num) = p_cent + dot_product(g_p_vec,Del_x_vec)
	d_ext_array(num) = d_cent + ((g_d_vec(1)*Del_x_vec(1))+(g_d_vec(2)*Del_x_vec(2))+(g_d_vec(3)*Del_x_vec(3)))
	ux_ext_array(num) = ux_cent + ((g_ux_vec(1)*Del_x_vec(1))+(g_ux_vec(2)*Del_x_vec(2))+(g_ux_vec(3)*Del_x_vec(3)))
	uy_ext_array(num) = uy_cent + ((g_uy_vec(1)*Del_x_vec(1))+(g_uy_vec(2)*Del_x_vec(2))+(g_uy_vec(3)*Del_x_vec(3)))
	uz_ext_array(num) = uz_cent + ((g_uz_vec(1)*Del_x_vec(1))+(g_uz_vec(2)*Del_x_vec(2))+(g_uz_vec(3)*Del_x_vec(3)))
	p_ext_array(num) = p_cent + ((g_p_vec(1)*Del_x_vec(1))+(g_p_vec(2)*Del_x_vec(2))+(g_p_vec(3)*Del_x_vec(3)))
!-------------------------------------------------------------------------------

!+ve x-face dire! phi_ngb,min = min_j[phi_j]ction neighbour
	idum = i+1
	jdum = j
	kdum = k

!Set position in array number:
	num=2

!Fill variable array entries for this neighbour:
! This is equivalent to entering the first value (j=1) into the phi_j, so that the min and max can be taken later.
	d_array(num) = cells(idum,jdum,kdum,rhon)
	ux_array(num) = cells(idum,jdum,kdum,uxn)
	uy_array(num) = cells(idum,jdum,kdum,uyn)
	uz_array(num) = cells(idum,jdum,kdum,uzn)
	p_array(num) = cells(idum,jdum,kdum,pn)


! *** EXTRAPOLATION ***:
!Calculate cell centre to cell face distances:
	delx = (cells(i,j,k,xn) - cells(idum,jdum,kdum,xn))/2.d0
	dely = (cells(i,j,k,yn) - cells(idum,jdum,kdum,yn))/2.d0	
	delz = (cells(i,j,k,zn) - cells(idum,jdum,kdum,zn))/2.d0	

!Fill temporary cell centre to face centre vector:
	Del_x_vec = (/delx,dely,delz/)

!Make extrapolations to cell face:
! " phi+grad_dum-vec·Del_x_j-vec "
!	d_ext_array(num) = d_cent + dot_product(g_d_vec,Del_x_vec)
!	ux_ext_array(num) = ux_cent + dot_product(g_ux_vec,Del_x_vec)
!	uy_ext_array(num) = uy_cent + dot_product(g_uy_vec,Del_x_vec)
!	uz_ext_array(num) = uz_cent + dot_product(g_uz_vec,Del_x_vec)
!	p_ext_array(num) = p_cent + dot_product(g_p_vec,Del_x_vec)
	d_ext_array(num) = d_cent + ((g_d_vec(1)*Del_x_vec(1))+(g_d_vec(2)*Del_x_vec(2))+(g_d_vec(3)*Del_x_vec(3)))
	ux_ext_array(num) = ux_cent + ((g_ux_vec(1)*Del_x_vec(1))+(g_ux_vec(2)*Del_x_vec(2))+(g_ux_vec(3)*Del_x_vec(3)))
	uy_ext_array(num) = uy_cent + ((g_uy_vec(1)*Del_x_vec(1))+(g_uy_vec(2)*Del_x_vec(2))+(g_uy_vec(3)*Del_x_vec(3)))
	uz_ext_array(num) = uz_cent + ((g_uz_vec(1)*Del_x_vec(1))+(g_uz_vec(2)*Del_x_vec(2))+(g_uz_vec(3)*Del_x_vec(3)))
	p_ext_array(num) = p_cent + ((g_p_vec(1)*Del_x_vec(1))+(g_p_vec(2)*Del_x_vec(2))+(g_p_vec(3)*Del_x_vec(3)))

!-------------------------------------------------------------------------------
!#####			y-Axis Neighbours				########
!-------------------------------------------------------------------------------

!-ve x-face direction neighbour
	idum = i
	jdum = j-1
	kdum = k

!Set position in array number:
	num=3

!Fill variable array entries for this neighbour:
! This is equivalent to entering the first value (j=1) into the phi_j, so that the min and max can be taken later.
	d_array(num) = cells(idum,jdum,kdum,rhon)
	ux_array(num) = cells(idum,jdum,kdum,uxn)
	uy_array(num) = cells(idum,jdum,kdum,uyn)
	uz_array(num) = cells(idum,jdum,kdum,uzn)
	p_array(num) = cells(idum,jdum,kdum,pn)


! *** EXTRAPOLATION ***:
!Calculate cell centre to cell face distances:
	delx = (cells(i,j,k,xn) - cells(idum,jdum,kdum,xn))/2.d0
	dely = (cells(i,j,k,yn) - cells(idum,jdum,kdum,yn))/2.d0	
	delz = (cells(i,j,k,zn) - cells(idum,jdum,kdum,zn))/2.d0	

!Fill temporary cell centre to face centre vector:
	Del_x_vec = (/delx,dely,delz/)

!Make extrapolations to cell face:
! " phi+grad_dum-vec·Del_x_j-vec "
!	d_ext_array(num) = d_cent + dot_product(g_d_vec,Del_x_vec)
!	ux_ext_array(num) = ux_cent + dot_product(g_ux_vec,Del_x_vec)
!	uy_ext_array(num) = uy_cent + dot_product(g_uy_vec,Del_x_vec)
!	uz_ext_array(num) = uz_cent + dot_product(g_uz_vec,Del_x_vec)
!	p_ext_array(num) = p_cent + dot_product(g_p_vec,Del_x_vec)
	d_ext_array(num) = d_cent + ((g_d_vec(1)*Del_x_vec(1))+(g_d_vec(2)*Del_x_vec(2))+(g_d_vec(3)*Del_x_vec(3)))
	ux_ext_array(num) = ux_cent + ((g_ux_vec(1)*Del_x_vec(1))+(g_ux_vec(2)*Del_x_vec(2))+(g_ux_vec(3)*Del_x_vec(3)))
	uy_ext_array(num) = uy_cent + ((g_uy_vec(1)*Del_x_vec(1))+(g_uy_vec(2)*Del_x_vec(2))+(g_uy_vec(3)*Del_x_vec(3)))
	uz_ext_array(num) = uz_cent + ((g_uz_vec(1)*Del_x_vec(1))+(g_uz_vec(2)*Del_x_vec(2))+(g_uz_vec(3)*Del_x_vec(3)))
	p_ext_array(num) = p_cent + ((g_p_vec(1)*Del_x_vec(1))+(g_p_vec(2)*Del_x_vec(2))+(g_p_vec(3)*Del_x_vec(3)))

!-------------------------------------------------------------------------------

!+ve x-face dire! phi_ngb,min = min_j[phi_j]ction neighbour
	idum = i
	jdum = j+1
	kdum = k

!Set position in array number:
	num=4

!Fill variable array entries for this neighbour:
! This is equivalent to entering the first value (j=1) into the phi_j, so that the min and max can be taken later.
	d_array(num) = cells(idum,jdum,kdum,rhon)
	ux_array(num) = cells(idum,jdum,kdum,uxn)
	uy_array(num) = cells(idum,jdum,kdum,uyn)
	uz_array(num) = cells(idum,jdum,kdum,uzn)
	p_array(num) = cells(idum,jdum,kdum,pn)


! *** EXTRAPOLATION ***:
!Calculate cell centre to cell face distances:
	delx = (cells(i,j,k,xn) - cells(idum,jdum,kdum,xn))/2.d0
	dely = (cells(i,j,k,yn) - cells(idum,jdum,kdum,yn))/2.d0	
	delz = (cells(i,j,k,zn) - cells(idum,jdum,kdum,zn))/2.d0	

!Fill temporary cell centre to face centre vector:
	Del_x_vec = (/delx,dely,delz/)

!Make extrapolations to cell face:
! " phi+grad_dum-vec·Del_x_j-vec "
!	d_ext_array(num) = d_cent + dot_product(g_d_vec,Del_x_vec)
!	ux_ext_array(num) = ux_cent + dot_product(g_ux_vec,Del_x_vec)
!	uy_ext_array(num) = uy_cent + dot_product(g_uy_vec,Del_x_vec)
!	uz_ext_array(num) = uz_cent + dot_product(g_uz_vec,Del_x_vec)
!	p_ext_array(num) = p_cent + dot_product(g_p_vec,Del_x_vec)
	d_ext_array(num) = d_cent + ((g_d_vec(1)*Del_x_vec(1))+(g_d_vec(2)*Del_x_vec(2))+(g_d_vec(3)*Del_x_vec(3)))
	ux_ext_array(num) = ux_cent + ((g_ux_vec(1)*Del_x_vec(1))+(g_ux_vec(2)*Del_x_vec(2))+(g_ux_vec(3)*Del_x_vec(3)))
	uy_ext_array(num) = uy_cent + ((g_uy_vec(1)*Del_x_vec(1))+(g_uy_vec(2)*Del_x_vec(2))+(g_uy_vec(3)*Del_x_vec(3)))
	uz_ext_array(num) = uz_cent + ((g_uz_vec(1)*Del_x_vec(1))+(g_uz_vec(2)*Del_x_vec(2))+(g_uz_vec(3)*Del_x_vec(3)))
	p_ext_array(num) = p_cent + ((g_p_vec(1)*Del_x_vec(1))+(g_p_vec(2)*Del_x_vec(2))+(g_p_vec(3)*Del_x_vec(3)))

!-------------------------------------------------------------------------------
!#####			z-Axis Neighbours				########
!-------------------------------------------------------------------------------

!-ve x-face direction neighbour
	idum = i
	jdum = j
	kdum = k-1

!Set position in array number:
	num=5

!Fill variable array entries for this neighbour:
! This is equivalent to entering the first value (j=1) into the phi_j, so that the min and max can be taken later.
	d_array(num) = cells(idum,jdum,kdum,rhon)
	ux_array(num) = cells(idum,jdum,kdum,uxn)
	uy_array(num) = cells(idum,jdum,kdum,uyn)
	uz_array(num) = cells(idum,jdum,kdum,uzn)
	p_array(num) = cells(idum,jdum,kdum,pn)


! *** EXTRAPOLATION ***:
!Calculate cell centre to cell face distances:
	delx = (cells(i,j,k,xn) - cells(idum,jdum,kdum,xn))/2.d0
	dely = (cells(i,j,k,yn) - cells(idum,jdum,kdum,yn))/2.d0	
	delz = (cells(i,j,k,zn) - cells(idum,jdum,kdum,zn))/2.d0	

!Fill temporary cell centre to face centre vector:
	Del_x_vec = (/delx,dely,delz/)

!Make extrapolations to cell face:
! " phi+grad_dum-vec·Del_x_j-vec "
!	d_ext_array(num) = d_cent + dot_product(g_d_vec,Del_x_vec)
!	ux_ext_array(num) = ux_cent + dot_product(g_ux_vec,Del_x_vec)
!	uy_ext_array(num) = uy_cent + dot_product(g_uy_vec,Del_x_vec)
!	uz_ext_array(num) = uz_cent + dot_product(g_uz_vec,Del_x_vec)
!	p_ext_array(num) = p_cent + dot_product(g_p_vec,Del_x_vec)
	d_ext_array(num) = d_cent + ((g_d_vec(1)*Del_x_vec(1))+(g_d_vec(2)*Del_x_vec(2))+(g_d_vec(3)*Del_x_vec(3)))
	ux_ext_array(num) = ux_cent + ((g_ux_vec(1)*Del_x_vec(1))+(g_ux_vec(2)*Del_x_vec(2))+(g_ux_vec(3)*Del_x_vec(3)))
	uy_ext_array(num) = uy_cent + ((g_uy_vec(1)*Del_x_vec(1))+(g_uy_vec(2)*Del_x_vec(2))+(g_uy_vec(3)*Del_x_vec(3)))
	uz_ext_array(num) = uz_cent + ((g_uz_vec(1)*Del_x_vec(1))+(g_uz_vec(2)*Del_x_vec(2))+(g_uz_vec(3)*Del_x_vec(3)))
	p_ext_array(num) = p_cent + ((g_p_vec(1)*Del_x_vec(1))+(g_p_vec(2)*Del_x_vec(2))+(g_p_vec(3)*Del_x_vec(3)))

!-------------------------------------------------------------------------------

!+ve x-face dire! phi_ngb,min = min_j[phi_j]ction neighbour
	idum = i
	jdum = j
	kdum = k+1

!Set position in array number:
	num=6

!Fill variable array entries for this neighbour:
! This is equivalent to entering the first value (j=1) into the phi_j, so that the min and max can be taken later.
	d_array(num) = cells(idum,jdum,kdum,rhon)
	ux_array(num) = cells(idum,jdum,kdum,uxn)
	uy_array(num) = cells(idum,jdum,kdum,uyn)
	uz_array(num) = cells(idum,jdum,kdum,uzn)
	p_array(num) = cells(idum,jdum,kdum,pn)


! *** EXTRAPOLATION ***:
!Calculate cell centre to cell face distances:
	delx = (cells(i,j,k,xn) - cells(idum,jdum,kdum,xn))/2.d0
	dely = (cells(i,j,k,yn) - cells(idum,jdum,kdum,yn))/2.d0	
	delz = (cells(i,j,k,zn) - cells(idum,jdum,kdum,zn))/2.d0	

!Fill temporary cell centre to face centre vector:
	Del_x_vec = (/delx,dely,delz/)

!Make extrapolations to cell face:
! " phi+grad_dum-vec·Del_x_j-vec "
!	d_ext_array(num) = d_cent + dot_product(g_d_vec,Del_x_vec)
!	ux_ext_array(num) = ux_cent + dot_product(g_ux_vec,Del_x_vec)
!	uy_ext_array(num) = uy_cent + dot_product(g_uy_vec,Del_x_vec)
!	uz_ext_array(num) = uz_cent + dot_product(g_uz_vec,Del_x_vec)
!	p_ext_array(num) = p_cent + dot_product(g_p_vec,Del_x_vec)
	d_ext_array(num) = d_cent + ((g_d_vec(1)*Del_x_vec(1))+(g_d_vec(2)*Del_x_vec(2))+(g_d_vec(3)*Del_x_vec(3)))
	ux_ext_array(num) = ux_cent + ((g_ux_vec(1)*Del_x_vec(1))+(g_ux_vec(2)*Del_x_vec(2))+(g_ux_vec(3)*Del_x_vec(3)))
	uy_ext_array(num) = uy_cent + ((g_uy_vec(1)*Del_x_vec(1))+(g_uy_vec(2)*Del_x_vec(2))+(g_uy_vec(3)*Del_x_vec(3)))
	uz_ext_array(num) = uz_cent + ((g_uz_vec(1)*Del_x_vec(1))+(g_uz_vec(2)*Del_x_vec(2))+(g_uz_vec(3)*Del_x_vec(3)))
	p_ext_array(num) = p_cent + ((g_p_vec(1)*Del_x_vec(1))+(g_p_vec(2)*Del_x_vec(2))+(g_p_vec(3)*Del_x_vec(3)))


!-------------------------------------------------------------------------------
!#####			FINAL CALCULATION AND GRADIENT ADJUST		########
!-------------------------------------------------------------------------------
!Minimum values of phi neighbours:
! phi_ngb,min = min_j[phi_j]
		d_ngb_min = minval(d_array)
		ux_ngb_min = minval(ux_array)
		uy_ngb_min = minval(uy_array)
		uz_ngb_min = minval(uz_array)
		p_ngb_min = minval(p_array)

!Maximum values of phi neighbours:
! phi_ngb,max = max_j[phi_j]
		d_ngb_max = maxval(d_array)
		ux_ngb_max = maxval(ux_array)
		uy_ngb_max = maxval(uy_array)
		uz_ngb_max = maxval(uz_array)
		p_ngb_max = maxval(p_array)

!Minimum values of phi EXTRAPOLATED:
! phi_ext,MIN = MIN_j[phi+grad_dum-vec·Del_x_j-vec]
		d_ext_min = minval(d_ext_array)
		ux_ext_min = minval(ux_ext_array)
		uy_ext_min = minval(uy_ext_array)
		uz_ext_min = minval(uz_ext_array)
		p_ext_min = minval(p_ext_array)

!Maximum values of phi EXTRAPOLATED:
! phi_ext,max = max_j[phi+grad_dum-vec·Del_x_j-vec]
		d_ext_max = maxval(d_ext_array)
		ux_ext_max = maxval(ux_ext_array)
		uy_ext_max = maxval(uy_ext_array)
		uz_ext_max = maxval(uz_ext_array)
		p_ext_max = maxval(p_ext_array)

!Calculate ALPHA for each variable and adjust gradient:
! alpha = minval (1,beta*minval[(phi_ngb,max - phi)/(phi_ext,max - phi), (phi - phi_ngb,min)/(phi - phi_ext,min)])
! Grad-vec = alpha* grad_dum-vec

!Density:
 alpha = max(0.d0,min(1.d0,beta*min((d_ngb_max - d_cent)/(d_ext_max - d_cent),(d_cent - d_ngb_min)/(d_cent - d_ext_min))))

 cells(i,j,k,gdxn) = alpha*g_d_vec(1)
 cells(i,j,k,gdyn) = alpha*g_d_vec(2)
 cells(i,j,k,gdzn) = alpha*g_d_vec(3)

!ux:
 alpha = max(0.d0,min(1.d0,beta*min((ux_ngb_max - ux_cent)/(ux_ext_max - ux_cent),(ux_cent - ux_ngb_min)/(ux_cent - ux_ext_min))))

 cells(i,j,k,guxxn) = alpha*g_ux_vec(1)
 cells(i,j,k,guxyn) = alpha*g_ux_vec(2)
 cells(i,j,k,guxzn) = alpha*g_ux_vec(3)


!uy:
 alpha = max(0.d0,min(1.d0,beta*min((uy_ngb_max - uy_cent)/(uy_ext_max - uy_cent),(uy_cent - uy_ngb_min)/(uy_cent - uy_ext_min))))

 cells(i,j,k,guyxn) = alpha*g_uy_vec(1)
 cells(i,j,k,guyyn) = alpha*g_uy_vec(2)
 cells(i,j,k,guyzn) = alpha*g_uy_vec(3)


!uz:
 alpha = max(0.d0,min(1.d0,beta*min((uz_ngb_max - uz_cent)/(uz_ext_max - uz_cent),(uz_cent - uz_ngb_min)/(uz_cent - uz_ext_min))))

 cells(i,j,k,guzxn) = alpha*g_uz_vec(1)
 cells(i,j,k,guzyn) = alpha*g_uz_vec(2)
 cells(i,j,k,guzzn) = alpha*g_uz_vec(3)


!p:
 alpha = max(0.d0,min(1.d0,beta*min((p_ngb_max - p_cent)/(p_ext_max - p_cent),(p_cent - p_ngb_min)/(p_cent - p_ext_min))))

 cells(i,j,k,gpxn) = alpha*g_p_vec(1)
 cells(i,j,k,gpyn) = alpha*g_p_vec(2)
 cells(i,j,k,gpzn) = alpha*g_p_vec(3)

end subroutine slope_limiter
