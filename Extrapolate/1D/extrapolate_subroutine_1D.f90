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
!* @file extrapolate_subroutine_1D.f90
!*
!* @Extrapolates primitive variables of cells from centre to cell edge for second
!* order functionality. 1D Version: Fortran 90 version.
!*
!* @author Andrew Hannington (ath4@st-andrews.ac.uk)
!********************************************************************************

!===============================================================================
!# Subroutine to extrapolate density, velocity components, and pressure from 
!# gradients calculated in "get gradients subroutine". Extrapolates the variable
!# to cell interface (face of cell with neighbour) and moves the variable forward
!# by half a timestep. This follows the equation:
!#
!# W_t + A(W)*W_x + B(W)*W_y + C(W)*W_z = 0
!#
!# Where:
!#
!# W=(/rho,ux,uy,uz,p/)
!#
!# W_t = Del/del(t) [W] etc.
!#
!# A(w)=(/(ux,rho,       0, 0,   0),
!#	  (0,  ux,       0, 0,1/rho),
!#	  (0,  0,        ux,0,   0),
!#	  (0,  0,        0, ux,  0),
!#	  (0,rho*(cs**2),0, 0,   ux)/)
!# B(w)=(/(uy,rho,  0,         0,   0),
!#	  (0,  uy,  0,         0,   0),
!#	  (0,  0,   uy        ,0,1/rho),
!#	  (0,  0,   0,         uy,  0),
!#	  (0,  0, rho*(cs**2), 0,   uy)/)
!# C(w)=(/(uz,rho,0,    0,         0),
!#	  (0,  uz,0,    0,         0),
!#	  (0,  0, uz,   0,         0),
!#	  (0,  0, 0,    uz,      1/rho),
!#	  (0,  0, 0,rho*(cs**2),   uz)/)
!#
!# Which leads to the half timestep extrapolation of:
!#
!# Delta[W] = -1.d0*(A(W)*W_x + B(W)*W_y + C(W)*W_z)*(Delta[t])/2.d0
!#
!===============================================================================
subroutine extrapolate (xyz,i,j,k,d1,u1vec,p1,d2,u2vec,p2,&
			& dL,uLvec,pL,dR,uRvec,pR,time)
	use constants
	implicit none

		Real*8, intent(in) :: time
		character(len=1), intent(in) :: xyz				!Axis indicator
		Integer, intent (in) :: i,j,k					!Indices of cell being analysed
		Real*8,intent(in):: d1,p1,d2,p2					!Density and pressure values in, for left and right
		Real*8, dimension(3),intent(in) :: u1vec,u2vec			!Velocity vector in, for left and right 
		Real*8,intent(out):: dL,pL,dR,pR				!extrapolated pressure and density at cell face, left and right
		Real*8, dimension(3),intent(out) :: uLvec,uRvec			!Extrapolated velocity at cell face, left and right
		Integer :: iL,iR,jL,jR,kL,kR					!Left and right values of i,j,k


!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)			

	If (xyz=="x") then
!===============================================================================
!#####				X-Axis Analysis				   #####
!===============================================================================
	! Set dummy left and right indices
	!	(x-axis):
		iL=i
		iR=i+1
		jL=j
		jR=j
		kL=k
		kR=k

		call extrapolate_calc (xyz,iL,iR,jL,jR,kL,kR,&
			& d1,u1vec,p1,d2,u2vec,p2,&
			& dL,uLvec,pL,dR,uRvec,pR, &
			& time)

	else if (xyz=="y") then
!===============================================================================
!#####				Y-axis Analysis				   #####
!===============================================================================
!!	! Set dummy left and right indices
!!	!	(y-axis):
!		iL=i
!		iR=i
!		jL=j
!		jR=j+1
!		kL=k
!		kR=k
!
!		call extrapolate_calc (xyz,iL,iR,jL,jR,kL,kR,&
!			& d1,u1vec,p1,d2,u2vec,p2,&
!			& dL,uLvec,pL,dR,uRvec,pR, &
!			& time)

	else if (xyz=="z") then
!===============================================================================
!#####				Z-Axis Analysis				   #####
!===============================================================================
	!! Set dummy left and right indices
	!!	(z-axis):
	!	iL=i
	!	iR=i
	!	jL=j
	!	jR=j
	!	kL=k
	!	kR=k+1
	!
!		call extrapolate_calc (xyz,iL,iR,jL,jR,kL,kR,&
!			& d1,u1vec,p1,d2,u2vec,p2,&
!			& dL,uLvec,pL,dR,uRvec,pR, &
!			& time)
		
	else
		Print*,"Axis indicator failure [@Extrapolate]! Program will Abort!"
		STOP
	endif


end subroutine extrapolate
!===============================================================================
!# Generalised subroutine to perform the extrapolation calculation. Only performs
!# spatial extrapolation on axis being analysed, and performs total temporal
!# extrapolation from above (main header) formulae.
!===============================================================================
subroutine extrapolate_calc (xyz,iL,iR,jL,jR,kL,kR,&
			& d1,u1vec,p1,d2,u2vec,p2,&
			& dLsol,uLvec,pLsol,dRsol,uRvec,pRsol, &
			& time)
	use constants
	implicit none
	
		Real*8, intent(in) :: time					!Current code time
		INTEGER,intent(in) :: iL,iR,jL,jR,kL,kR				!Left and right indices
		Character(len=1), intent(in) :: xyz

		Real*8,intent(in):: d1,p1,d2,p2					!Density and pressure values in, for left and right
		Real*8, dimension(3),intent(in) :: u1vec,u2vec			!Velocity vector in, for left and right 
		Real*8,intent(out):: dLsol,pLsol,dRsol,pRsol			!extrapolated pressure and density at cell face, left and right
		Real*8, dimension(3),intent(out) :: uLvec,uRvec			!Extrapolated velocity at cell face, left and right

		Real*8 :: xL,xR,yL,yR,zL,zR					!cell position for left and right cell, for given axis direction
		Real*8 :: dL,uxL,uyL,uzL,pL					!Left cell quantities. This cell will be the cell under analysis
		Real*8 :: dR,uxR,uyR,uzR,pR					!Right cell quantities.

		Real*8 :: grad_d_x_L,grad_d_y_L,grad_d_z_L, &			!Left Density gradients on each axis
			& grad_d_x_R,grad_d_y_R,grad_d_z_R			!Right Density gradients on each axis

		Real*8 :: grad_ux_x_L,grad_uy_x_L,grad_uz_x_L, &		!Left velocity gradients for each component on each axis
			& grad_ux_y_L,grad_uy_y_L,grad_uz_y_L, &
			& grad_ux_z_L,grad_uy_z_L,grad_uz_z_L
		Real*8 :: grad_ux_x_R,grad_uy_x_R,grad_uz_x_R, &		!Right velocity gradients for each component on each axis
			& grad_ux_y_R,grad_uy_y_R,grad_uz_y_R, &
			& grad_ux_z_R,grad_uy_z_R,grad_uz_z_R

		Real*8 :: grad_p_x_L,grad_p_y_L,grad_p_z_L, &			!Left pressure gradients on each axis
			& grad_p_x_R,grad_p_y_R,grad_p_z_R			!Right pressure gradients on each axis

		Real*8 :: d_extrap_L,ux_extrap_L,uy_extrap_L,uz_extrap_L,p_extrap_L, &
			& d_extrap_R,ux_extrap_R,uy_extrap_R,uz_extrap_R,p_extrap_R	!Exratoplated variables, left and right

		Real*8 :: Delx,Dely,Delz					!Cell-to-cell distance	on each axis
	
		Real*8 :: cs_L, cs_R						!Sound speed in left and rigt cell cs=Sqrt(gamma*(p/rho))		

		Real*8, dimension(5)	:: W_L, W_x_L, W_y_L, W_z_L, &		!Left and right W vectors, and spatial derivatives
					&  W_R, W_x_R, W_y_R, W_z_R

		Real*8, dimension (5,5)	:: AMat_L, BMat_L, CMat_L, &		!Left and right A,B,C matrices
					 & AMat_R, BMat_R, CMat_R
		Integer :: h							!Counter for do loop
		Character(len=4) :: error_loc					!Error location identifier string		




!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)

!-------------------------------------------------------------------------------
!# Set variable values from cells being analysed:
!-------------------------------------------------------------------------------
!	Set variables from cells to be used:
!		Left values:
!		x=:
	xL=cells(iL,jL,kL,xn)
!		y=:
	yL=cells(iL,jL,kL,yn)
!		z=:
	zL=cells(iL,jL,kL,zn)

!	Density:
	dL=d1

!	Pressure:
	pL=p1

!	Velocity components:
	uxL = u1vec(1)
	uyL = u1vec(2)
	uzL = u1vec(3)

!Left Density gradients on each axis:
	grad_d_x_L=cells(iL,jL,kL,gdxn)
	grad_d_y_L=cells(iL,jL,kL,gdyn)
	grad_d_z_L=cells(iL,jL,kL,gdzn)			

!Left velocity gradients for each component on each axis:
	grad_ux_x_L=cells(iL,jL,kL,guxxn)
	grad_uy_x_L=cells(iL,jL,kL,guyxn)
	grad_uz_x_L=cells(iL,jL,kL,guzxn)
	
	grad_ux_y_L=cells(iL,jL,kL,guxyn)
	grad_uy_y_L=cells(iL,jL,kL,guyyn)
	grad_uz_y_L=cells(iL,jL,kL,guzyn)

	grad_ux_z_L=cells(iL,jL,kL,guxzn)
	grad_uy_z_L=cells(iL,jL,kL,guyzn)
	grad_uz_z_L=cells(iL,jL,kL,guzzn)

!Left pressure gradients on each axis:
	grad_p_x_L=cells(iL,jL,kL,gpxn)
	grad_p_y_L=cells(iL,jL,kL,gpyn)
	grad_p_z_L=cells(iL,jL,kL,gpzn)

!-------------------------------------------------------------------------------
!	Set variables from cells to be used:
!		Right values:
!		x=:
	xR=cells(iR,jR,kR,xn)
!		y=:
	yR=cells(iR,jR,kR,yn)
!		z=:
	zR=cells(iR,jR,kR,zn)

!	Density:
	dR=d2

!	Pressure:
	pR=p2

!	Velocity components:
	uxR = u2vec(1)
	uyR = u2vec(2)
	uzR = u2vec(3)

!Left Density gradients on each axis:
	grad_d_x_R=cells(iR,jR,kR,gdxn)
	grad_d_y_R=cells(iR,jR,kR,gdyn)
	grad_d_z_R=cells(iR,jR,kR,gdzn)			

!Left velocity gradients for each component on each axis:
	grad_ux_x_R=cells(iR,jR,kR,guxxn)
	grad_uy_x_R=cells(iR,jR,kR,guyxn)
	grad_uz_x_R=cells(iR,jR,kR,guzxn)
	
	grad_ux_y_R=cells(iR,jR,kR,guxyn)
	grad_uy_y_R=cells(iR,jR,kR,guyyn)
	grad_uz_y_R=cells(iR,jR,kR,guzyn)

	grad_ux_z_R=cells(iR,jR,kR,guxzn)
	grad_uy_z_R=cells(iR,jR,kR,guyzn)
	grad_uz_z_R=cells(iR,jR,kR,guzzn)

!Left pressure gradients on each axis:
	grad_p_x_R=cells(iR,jR,kR,gpxn)
	grad_p_y_R=cells(iR,jR,kR,gpyn)
	grad_p_z_R=cells(iR,jR,kR,gpzn)

!-------------------------------------------------------------------------------
!# Sets extrapolation direction for relevant axis being analysed. Thus, centre
!# to relevant cell face on given axis, distance is calculated for extrapolation.
!-------------------------------------------------------------------------------
	If (xyz=="x") then
		Delx=DABS((xR-xL)/2.d0)
		Dely=0.d0
		Delz=0.d0
	else if (xyz=="y") then
		Delx=0.d0
		Dely=DABS((yR-yL)/2.d0)
		Delz=0.d0
	else if (xyz=="z") then
		Delx=0.d0
		Dely=0.d0
		Delz=DABS((zR-zL)/2.d0)
	else
		Print*,"Axis indicator failure [@Extrapolate_calc]! Program will Abort!"
		STOP
	endif

	If (xyz=="x") then
!-------------------------------------------------------------------------------
!# 			x-axis Extrapolation
!-------------------------------------------------------------------------------
!	Extrapolate values of pressure, density, and velocity:
!		Left Cell:
		d_extrap_L = dL + Delx*grad_d_x_L
		ux_extrap_L = uxL + Delx*grad_ux_x_L
		uy_extrap_L = uyL + Delx*grad_uy_x_L
		uz_extrap_L = uzL + Delx*grad_uz_x_L
		p_extrap_L = pL + Delx*grad_p_x_L
!		Right Cell:
		d_extrap_R = dR - Delx*grad_d_x_R
		ux_extrap_R = uxR - Delx*grad_ux_x_R
		uy_extrap_R = uyR - Delx*grad_uy_x_R
		uz_extrap_R = uzR - Delx*grad_uz_x_R
		p_extrap_R = pR - Delx*grad_p_x_R
	else if (xyz=="y") then
!-------------------------------------------------------------------------------
!# 			y-axis Extrapolation
!-------------------------------------------------------------------------------
!	Extrapolate values of pressure, density, and velocity:
!		Left Cell:
		d_extrap_L = dL + Dely*grad_d_y_L
		ux_extrap_L = uxL + Dely*grad_ux_y_L
		uy_extrap_L = uyL + Dely*grad_uy_y_L
		uz_extrap_L = uzL + Dely*grad_uz_y_L
		p_extrap_L = pL + Dely*grad_p_y_L
!		Right Cell:
		d_extrap_R = dR - Dely*grad_d_y_R
		ux_extrap_R = uxR - Dely*grad_ux_y_R
		uy_extrap_R = uyR - Dely*grad_uy_y_R
		uz_extrap_R = uzR - Dely*grad_uz_y_R
		p_extrap_R = pR - Dely*grad_p_y_R
	else if (xyz=="z") then
!-------------------------------------------------------------------------------
!# 			z-axis Extrapolation
!-------------------------------------------------------------------------------
!	Extrapolate values of pressure, density, and velocity:
!		Left Cell:
		d_extrap_L = dL + Delz*grad_d_z_L
		ux_extrap_L = uxL + Delz*grad_ux_z_L
		uy_extrap_L = uyL + Delz*grad_uy_z_L
		uz_extrap_L = uzL + Delz*grad_uz_z_L
		p_extrap_L = pL + Delz*grad_p_z_L
!		Right Cell:
		d_extrap_R = dR - Delz*grad_d_z_R
		ux_extrap_R = uxR - Delz*grad_ux_z_R
		uy_extrap_R = uyR - Delz*grad_uy_z_R
		uz_extrap_R = uzR - Delz*grad_uz_z_R
		p_extrap_R = pR - Delz*grad_p_z_R
	else
		Print*,"Axis indicator failure [@Extrapolate_calc]! Program will Abort!"
		STOP
	endif

!	Check for unphysical values:
	error_loc="extr"

!	call error_sub(0,error_loc,NcXr,NcYr,NcZr,gam,1.d0,1.d0,cellVol,p_extrap_L,p_extrap_R,d_extrap_L,d_extrap_R,iL,jL,kL,time)

!-------------------------------------------------------------------------------
!# 			Setup A, B, C matrices, and W vectors:
!-------------------------------------------------------------------------------
! LEFT and RIGHT W vectors
	W_L(1:5) = (/d_extrap_L,ux_extrap_L,uy_extrap_L,uz_extrap_L,p_extrap_L/)
	W_R(1:5) = (/d_extrap_R,ux_extrap_R,uy_extrap_R,uz_extrap_R,p_extrap_R/)
! LEFT and RIGHT W vector spatial derivatives:
	!Left cell:
	W_x_L(1:5) = (/grad_d_x_L,grad_ux_x_L,grad_uy_x_L,grad_uz_x_L,grad_p_x_L/)
	W_y_L(1:5) = (/grad_d_y_L,grad_ux_y_L,grad_uy_y_L,grad_uz_y_L,grad_p_y_L/)
	W_z_L(1:5) = (/grad_d_z_L,grad_ux_z_L,grad_uy_z_L,grad_uz_z_L,grad_p_z_L/)

	!Right cell:
	W_x_R(1:5) = (/grad_d_x_R,grad_ux_x_R,grad_uy_x_R,grad_uz_x_R,grad_p_x_R/)
	W_y_R(1:5) = (/grad_d_y_R,grad_ux_y_R,grad_uy_y_R,grad_uz_y_R,grad_p_y_R/)
	W_z_R(1:5) = (/grad_d_z_R,grad_ux_z_R,grad_uy_z_R,grad_uz_z_R,grad_p_z_R/)

! Left A matrix:
	AMat_L(1:5,1:5) = transpose(reshape( &						!Transpose and reshape takes the list as written, turns it into a matrix of the shape given by shape(AMat_L),
			  & (/uxL,dL,0.d0,0.d0,0.d0, &					! and transposes so the matrix output is of the format in the comments above and as written roughly here.
			  & 0.d0,uxL,0.d0,0.d0,(1.d0/dL), &
			  & 0.d0,0.d0,uxL,0.d0,0.d0, &
			  & 0.d0,0.d0,0.d0,uxL,0.d0, &
			  & 0.d0,gam*pL,0.d0,0.d0,uxL/),shape(AMat_L)))
! RIGHT A matrix:
	AMat_R(1:5,1:5) = transpose(reshape( &
			  & (/uxR,dR,0.d0,0.d0,0.d0, &
			  & 0.d0,uxR,0.d0,0.d0,(1.d0/dR), &
			  & 0.d0,0.d0,uxR,0.d0,0.d0, &
			  & 0.d0,0.d0,0.d0,uxR,0.d0, &
			  & 0.d0,gam*pR,0.d0,0.d0,uxR/),shape(AMat_R)))
! Left B matrix:
	BMat_L(1:5,1:5) = transpose(reshape( &
			  & (/uyL,dL,0.d0,0.d0,0.d0, &
			  & 0.d0,uyL,0.d0,0.d0,0.d0, &
			  & 0.d0,0.d0,uyL,0.d0,(1.d0/dL), &
			  & 0.d0,0.d0,0.d0,uyL,0.d0, &
			  & 0.d0,0.d0,gam*pL,0.d0,uyL/),shape(BMat_L)))
! RIGHT B matrix:
	BMat_R(1:5,1:5) = transpose(reshape( &
			  & (/uyR,dR,0.d0,0.d0,0.d0, &
			  & 0.d0,uyR,0.d0,0.d0,0.d0, &
			  & 0.d0,0.d0,uyR,0.d0,(1.d0/dR), &
			  & 0.d0,0.d0,0.d0,uyR,0.d0, &
			  & 0.d0,0.d0,gam*pR,0.d0,uyR/),shape(BMat_R)))
! Left C matrix:
	CMat_L(1:5,1:5) = transpose(reshape( &
			  & (/uzL,dL,0.d0,0.d0,0.d0, &
			  & 0.d0,uzL,0.d0,0.d0,0.d0, &
			  & 0.d0,0.d0,uzL,0.d0,0.d0, &
			  & 0.d0,0.d0,0.d0,uzL,(1.d0/dL), &
			  & 0.d0,0.d0,0.d0,gam*pL,uzL/),shape(CMat_L)))
! Left C matrix:
	CMat_R(1:5,1:5) = transpose(reshape( &
			  & (/uzR,dR,0.d0,0.d0,0.d0, &
			  & 0.d0,uzR,0.d0,0.d0,0.d0, &
			  & 0.d0,0.d0,uzR,0.d0,0.d0, &
			  & 0.d0,0.d0,0.d0,uzR,(1.d0/dR), &
			  & 0.d0,0.d0,0.d0,gam*pR,uzR/),shape(CMat_R)))

!-------------------------------------------------------------------------------
!# 			Predict forward by one half-timestep:
!-------------------------------------------------------------------------------
! LEFT predictions of quantities by half-timestep:

	h=1
	Do while (h<=5)								!Scrolls over each row updating each seperately. This allows for the dot_product function to be used.
		W_L(h) = W_L(h) - delta_t*0.5d0*( &
			& (AMat_L(h,1)*W_x_L(1)+AMat_L(h,2)*W_x_L(2) + &
			& AMat_L(h,3)*W_x_L(3)+AMat_L(h,4)*W_x_L(4) + &
			& AMat_L(h,5)*W_x_L(5)) + &
			& (BMat_L(h,1)*W_y_L(1)+BMat_L(h,2)*W_y_L(2) + &
			& BMat_L(h,3)*W_y_L(3)+BMat_L(h,4)*W_y_L(4) + &
			& BMat_L(h,3)*W_y_L(3)) + &
			& (CMat_L(h,1)*W_z_L(1)+CMat_L(h,2)*W_z_L(2) + &
			& CMat_L(h,3)*W_z_L(3)+CMat_L(h,4)*W_z_L(4) + &
			& CMat_L(h,5)*W_z_L(5)))

!		W_L(h) = W_L(h) - delta_t*0.5d0*(&
!			& dot_product(AMat_L(h,1:5),W_x_L(1:5)) + &
!			& dot_product(BMat_L(h,1:5),W_y_L(1:5)) + &
!			& dot_product(CMat_L(h,1:5),W_z_L(1:5)))
		h=h+1
	enddo

! RIGHT predictions of quantities by half-timestep:
	h=1
	Do while (h<=5)
		W_R(h) = W_R(h) - delta_t*0.5d0*( &
			& (AMat_R(h,1)*W_x_R(1)+AMat_R(h,2)*W_x_R(2) + &
			& AMat_R(h,3)*W_x_R(3)+AMat_R(h,4)*W_x_R(4) + &
			& AMat_R(h,5)*W_x_R(5)) + &
			& (BMat_R(h,1)*W_y_R(1)+BMat_R(h,2)*W_y_R(2) + &
			& BMat_R(h,3)*W_y_R(3)+BMat_R(h,4)*W_y_R(4) + &
			& BMat_R(h,5)*W_y_R(5)) + &
			& (CMat_R(h,1)*W_z_R(1)+CMat_R(h,2)*W_z_R(2) + &
			& CMat_R(h,3)*W_z_R(3)+CMat_R(h,4)*W_z_R(4) + &
			& CMat_R(h,5)*W_z_R(5)))

!		W_R(h) = W_R(h) - delta_t*0.5d0*(&
!			& dot_product(AMat_R(h,1:5),W_x_R(1:5)) + &
!			& dot_product(BMat_R(h,1:5),W_y_R(1:5)) + &
!			& dot_product(CMat_R(h,1:5),W_z_R(1:5)))
		h=h+1
	enddo	

!-------------------------------------------------------------------------------
!#				Update cells:
!-------------------------------------------------------------------------------

!		Left values:
!			Rho=:
		dLsol = W_L(1)
!			Velocity:
		uLvec(1:3) = (/W_L(2),W_L(3),W_L(4)/)
!			Pressure (p)=:
		pLsol = W_L(5)
!
!		Right values:
!			Rho=:
		dRsol = W_R(1)
!			Velocity:
		uRvec(1:3) = (/W_R(2),W_R(3),W_R(4)/)
!			Pressure (p)=:
		pRsol = W_R(5)

!	Check for unphysical values:
	error_loc="extr"
	call error_sub(0,error_loc,NcXr,NcYr,NcZr,gam,1.d0,1.d0,cellVol,W_L(5),W_R(5),W_L(1),W_R(1),iL,jL,kL,time)

end subroutine extrapolate_calc
